from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Set, Optional

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
    """Group SNP columns by the GWAS contrast they came from.

    - snp_columns = all merged columns that start with snp_prefix (default: 'rs')
    - comparison_group = sorted(Case, Control) joined by '_vs_'
    - for each comparison group, take meta SNP rsIDs and match merged SNP columns that startwith rsid + '_'
    """
    snp_columns = [col for col in merged.columns if str(col).startswith(snp_prefix)]

    meta = snp_metadata.copy()
    if normalize_case_control:
        meta["Case"] = meta["Case"].map(normalize_subtype)
        meta["Control"] = meta["Control"].map(normalize_subtype)

    meta["comparison_group"] = meta.apply(
        lambda row: "_vs_".join(sorted([str(row["Case"]), str(row["Control"])])),
        axis=1,
    )

    groups: Dict[str, List[str]] = {}
    for group_name in meta["comparison_group"].unique():
        rs_ids_for_group = meta.loc[meta["comparison_group"] == group_name, "SNP"].astype(str).tolist()
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
    """Run OLS for each SNP × connectivity pair.

    Model:
        connectivity ~ snp * C(Subtype, Treatment(reference)) + age + sex + C(scanner) + C(scan_type)

    Outputs per row:
      - R_squared
      - p-values / coefficients for each non-reference subtype vs reference (interaction terms)
      - p-values / coefficients for pairwise contrasts among non-reference subtypes (via t_test)
    """
    connectivity_columns = [c for c in merged.columns if connectivity_col_contains in str(c)]
    if not connectivity_columns:
        raise ValueError(
            f"No connectivity columns found matching pattern '{connectivity_col_contains}'. "
            f"Tip: use --connectivity-col-contains to change the pattern."
        )

    if subtype_categories is None:
        # infer from data (ordered categorical if set)
        subtype_categories = [str(x) for x in pd.Series(merged["Subtype"].dropna().unique()).tolist()]

    # non-reference subtypes
    nonref = [s for s in subtype_categories if s != reference_subtype]
    if len(nonref) < 1:
        raise ValueError("Need at least one non-reference subtype to fit interaction models.")

    # flatten all SNPs
    all_gwas_snps: Set[str] = {snp for snps in gwas_snp_groups.values() for snp in snps}

    formula_tmpl = "{conn} ~ {snp} * C(Subtype, Treatment('{ref}')) + age + sex + C(scanner) + C(scan_type)"
    results: List[dict] = []

    total = len(all_gwas_snps) * len(connectivity_columns)
    iterator = tqdm(all_gwas_snps, desc="SNPs", unit="snp") if progress else all_gwas_snps
    processed = 0

    for snp_col in iterator:
        for conn_col in connectivity_columns:
            processed += 1
            cols = [conn_col, snp_col, "Subtype", "age", "sex", "scanner", "scan_type"]
            try:
                subset = merged[cols].dropna()
            except KeyError:
                continue

            if len(subset) < min_obs_for_model:
                continue
            if subset[snp_col].std() == 0 or subset[conn_col].std() == 0:
                continue
            if subset["Subtype"].nunique(dropna=True) < 2:
                continue

            row = {"SNP_Name": snp_col, "Connectivity_Name": conn_col, "R_squared": np.nan}

            gwas_comp = [name for name, snps in gwas_snp_groups.items() if snp_col in snps]
            row["GWAS_Comparison"] = gwas_comp[0] if gwas_comp else "N/A"

            try:
                model = smf.ols(
                    formula_tmpl.format(conn=conn_col, snp=snp_col, ref=reference_subtype),
                    data=subset,
                ).fit()

                param_names = model.params.index.tolist()
                row["R_squared"] = float(model.rsquared)

                # non-reference vs reference interaction terms
                for subtype in nonref:
                    pname = f"{snp_col}:C(Subtype, Treatment('{reference_subtype}'))[T.{subtype}]"
                    if pname in param_names:
                        row[f"P_{subtype}_vs_{reference_subtype}"] = float(model.pvalues[pname])
                        row[f"Coeff_{subtype}_vs_{reference_subtype}"] = float(model.params[pname])

                # pairwise non-reference contrasts
                for g1, g2 in combinations(nonref, 2):
                    p1 = f"{snp_col}:C(Subtype, Treatment('{reference_subtype}'))[T.{g1}]"
                    p2 = f"{snp_col}:C(Subtype, Treatment('{reference_subtype}'))[T.{g2}]"
                    if p1 in param_names and p2 in param_names:
                        contrast = np.zeros(len(param_names))
                        contrast[param_names.index(p1)] = 1.0
                        contrast[param_names.index(p2)] = -1.0
                        t_res = model.t_test(contrast)
                        p_val = t_res.pvalue.item() if isinstance(t_res.pvalue, np.ndarray) else float(t_res.pvalue)
                        eff = t_res.effect.item() if isinstance(t_res.effect, np.ndarray) else float(t_res.effect)
                        row[f"P_{g1}_vs_{g2}"] = p_val
                        row[f"Coeff_{g1}_vs_{g2}"] = eff

                results.append(row)
            except Exception:
                pass

            if (not progress) and (processed % progress_every == 0):
                print(f"Processed {processed}/{total} combinations...")

    return pd.DataFrame(results)
