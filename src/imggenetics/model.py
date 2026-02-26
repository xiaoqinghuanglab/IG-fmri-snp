"""Core modeling routines for SNP × subtype interaction analysis.

The key analysis step fits an OLS model for every (SNP, connectivity-edge) pair.
We follow the notebook design:

  connectivity_edge ~ SNP * C(Subtype, Treatment(reference)) + covariates

The interaction terms encode subtype-specific SNP slopes relative to the
reference subtype (typically CN/Control). We then compute pairwise contrasts
between non-reference subtypes using model-based t-tests.

This module is written to be:
  - easy to read
  - robust to missing data (drops NAs per test)
  - safe for batch runs (skips combinations with insufficient variance)
"""

from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Optional, Set

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from tqdm import tqdm

from .utils import normalize_subtype


def build_gwas_snp_groups(
    merged: pd.DataFrame,
    snp_metadata: pd.DataFrame,
    snp_prefix: str = "rs",
    normalize_case_control: bool = True,
) -> Dict[str, List[str]]:
    """Group SNP columns by the GWAS contrast they originated from.

    Parameters
    ----------
    merged:
        The merged modeling table containing SNP columns.
    snp_metadata:
        A table with (at least) columns: Case, Control, SNP.
        Each row indicates a prioritized SNP from a given GWAS comparison.
    snp_prefix:
        Prefix used to identify SNP columns in `merged`.
    normalize_case_control:
        If True, normalize Case/Control labels in metadata using normalize_subtype.

    Returns
    -------
    dict
        Mapping: <comparison_group> -> list of SNP columns present in merged.

    Notes
    -----
    - comparison_group is built as: "<A>_vs_<B>" with A/B sorted alphabetically.
      This keeps naming stable regardless of direction in the metadata.
    - Genotype-matrix exports often encode SNP columns as "rs123_A" etc.
      We therefore match by prefix: rsID + "_".
    """

    # All SNP columns in the merged data.
    snp_columns = [col for col in merged.columns if str(col).startswith(snp_prefix)]

    meta = snp_metadata.copy()
    if normalize_case_control:
        meta["Case"] = meta["Case"].map(normalize_subtype)
        meta["Control"] = meta["Control"].map(normalize_subtype)

    # Stable comparison naming.
    meta["comparison_group"] = meta.apply(
        lambda row: "_vs_".join(sorted([str(row["Case"]), str(row["Control"])])), axis=1
    )

    groups: Dict[str, List[str]] = {}
    for group_name in meta["comparison_group"].unique():
        rs_ids_for_group = meta.loc[meta["comparison_group"] == group_name, "SNP"].astype(str).tolist()

        # Match SNP columns in the merged data. Many matrices include allele suffixes.
        actual_snp_cols_in_data = [
            col for col in snp_columns if any(str(col).startswith(rs_id + "_") for rs_id in rs_ids_for_group)
        ]
        groups[group_name] = actual_snp_cols_in_data

    return groups


def run_ols_interaction_grid(
    merged: pd.DataFrame,
    gwas_snp_groups: Dict[str, List[str]],
    connectivity_col_contains: str = "Connectivity",
    reference_subtype: str = "Control",
    subtype_categories: Optional[List[str]] = None,
    min_obs_for_model: int = 10,
    progress: bool = True,
    progress_every: int = 100000,
) -> pd.DataFrame:
    """Fit OLS models for all SNP × edge combinations.

    Model (formula)
    --------------
    connectivity_edge ~ snp * C(Subtype, Treatment(reference)) + age + sex + C(scanner) + C(scan_type)

    Output columns
    --------------
    Each output row corresponds to one (SNP, edge) test and includes:
      - SNP_Name, Connectivity_Name, GWAS_Comparison
      - R_squared
      - P_* and Coeff_* columns for:
          (non-reference subtype) vs (reference subtype)
      - P_* and Coeff_* for pairwise comparisons among non-reference subtypes
        computed via t_test (difference of interaction coefficients).

    Filtering / skipping
    --------------------
    We skip tests when:
      - too few complete observations
      - SNP or edge has zero variance
      - only one subtype present after NA-drop

    This keeps batch runs from failing on corner cases.
    """

    # Identify connectivity features.
    connectivity_columns = [c for c in merged.columns if connectivity_col_contains in str(c)]
    if not connectivity_columns:
        raise ValueError(
            f"No connectivity columns found matching pattern '{connectivity_col_contains}'. "
            f"Tip: use --connectivity-col-contains to change the pattern."
        )

    # Determine subtype categories (ordered categorical is best, but we can infer).
    if subtype_categories is None:
        subtype_categories = [str(x) for x in pd.Series(merged["Subtype"].dropna().unique()).tolist()]

    nonref = [s for s in subtype_categories if s != reference_subtype]
    if len(nonref) < 1:
        raise ValueError("Need at least one non-reference subtype to fit interaction models.")

    # Flatten all SNPs across groups (unique set).
    all_gwas_snps: Set[str] = {snp for snps in gwas_snp_groups.values() for snp in snps}

    # Covariates are included to match the original analysis intent.
    formula_tmpl = "{conn} ~ {snp} * C(Subtype, Treatment('{ref}')) + age + sex + C(scanner) + C(scan_type)"

    results: List[dict] = []

    total = len(all_gwas_snps) * len(connectivity_columns)
    iterator = tqdm(all_gwas_snps, desc="SNPs", unit="snp") if progress else all_gwas_snps
    processed = 0

    for snp_col in iterator:
        for conn_col in connectivity_columns:
            processed += 1

            # Work on complete-case rows for this specific test.
            cols = [conn_col, snp_col, "Subtype", "age", "sex", "scanner", "scan_type"]
            try:
                subset = merged[cols].dropna()
            except KeyError:
                # If any required covariate column is missing, skip.
                continue

            # Quick guards for degenerate cases.
            if len(subset) < min_obs_for_model:
                continue
            if subset[snp_col].std() == 0 or subset[conn_col].std() == 0:
                continue
            if subset["Subtype"].nunique(dropna=True) < 2:
                continue

            row = {"SNP_Name": snp_col, "Connectivity_Name": conn_col, "R_squared": np.nan}

            # Store which GWAS contrast produced this SNP (for downstream grouping).
            gwas_comp = [name for name, snps in gwas_snp_groups.items() if snp_col in snps]
            row["GWAS_Comparison"] = gwas_comp[0] if gwas_comp else "N/A"

            try:
                model = smf.ols(
                    formula_tmpl.format(conn=conn_col, snp=snp_col, ref=reference_subtype),
                    data=subset,
                ).fit()

                param_names = model.params.index.tolist()
                row["R_squared"] = float(model.rsquared)

                # 1) Non-reference vs reference interaction terms
                for subtype in nonref:
                    pname = f"{snp_col}:C(Subtype, Treatment('{reference_subtype}'))[T.{subtype}]"
                    if pname in param_names:
                        row[f"P_{subtype}_vs_{reference_subtype}"] = float(model.pvalues[pname])
                        row[f"Coeff_{subtype}_vs_{reference_subtype}"] = float(model.params[pname])

                # 2) Pairwise non-reference contrasts: (g1 - g2)
                for g1, g2 in combinations(nonref, 2):
                    p1 = f"{snp_col}:C(Subtype, Treatment('{reference_subtype}'))[T.{g1}]"
                    p2 = f"{snp_col}:C(Subtype, Treatment('{reference_subtype}'))[T.{g2}]"
                    if p1 in param_names and p2 in param_names:
                        contrast = np.zeros(len(param_names))
                        contrast[param_names.index(p1)] = 1.0
                        contrast[param_names.index(p2)] = -1.0

                        t_res = model.t_test(contrast)

                        # statsmodels sometimes returns scalars or 1x1 arrays.
                        p_val = t_res.pvalue.item() if isinstance(t_res.pvalue, np.ndarray) else float(t_res.pvalue)
                        eff = t_res.effect.item() if isinstance(t_res.effect, np.ndarray) else float(t_res.effect)

                        row[f"P_{g1}_vs_{g2}"] = p_val
                        row[f"Coeff_{g1}_vs_{g2}"] = eff

                results.append(row)

            except Exception:
                # We intentionally swallow errors here to keep long runs from failing.
                # If you want to debug, set progress=False and insert prints/raise.
                pass

            if (not progress) and (processed % progress_every == 0):
                print(f"Processed {processed}/{total} combinations...")

    return pd.DataFrame(results)
