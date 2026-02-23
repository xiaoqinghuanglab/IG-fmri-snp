from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from itertools import combinations
from typing import Dict, Iterable, List, Optional, Tuple, Set

from tqdm import tqdm



def build_gwas_snp_groups(
    merged: pd.DataFrame,
    snp_metadata: pd.DataFrame,
    snp_prefix: str = "rs",
) -> Dict[str, List[str]]:
    """
    Reproduce notebook logic:
      - snp_columns = all merged columns that start with 'rs'
      - create comparison_group = sorted(Case, Control) joined by '_vs_'
      - for each comparison group, take meta SNP rsIDs and match merged SNP columns that startwith rsid + '_'
    """
    snp_columns = [col for col in merged.columns if str(col).startswith(snp_prefix)]

    # compute comparison_group
    snp_metadata = snp_metadata.copy()
    snp_metadata["comparison_group"] = snp_metadata.apply(
        lambda row: "_vs_".join(sorted([str(row["Case"]), str(row["Control"])])),
        axis=1
    )

    groups: Dict[str, List[str]] = {}
    for group_name in snp_metadata["comparison_group"].unique():
        rs_ids_for_group = snp_metadata.loc[snp_metadata["comparison_group"] == group_name, "SNP"].astype(str).tolist()
        actual_snp_cols_in_data = [
            col for col in snp_columns
            if any(str(col).startswith(rs_id + "_") for rs_id in rs_ids_for_group)
        ]
        groups[group_name] = actual_snp_cols_in_data
    return groups


def run_ols_interaction_grid(
    merged: pd.DataFrame,
    gwas_snp_groups: Dict[str, List[str]],
    connectivity_col_contains: str = "Connectivity",
    min_obs_for_model: int = 10,
    progress: bool = True,
    progress_every: int = 100000,
) -> pd.DataFrame:
    """
    Run OLS for each SNP × connectivity pair, mirroring notebook logic.

    Model:
      connectivity ~ snp * C(Subtype, Treatment('Control')) + age + sex + C(scanner) + C(scan_type)

    Outputs:
      rows with p-values and coefficients for:
        - subtype vs Control (interaction terms)
        - pairwise subtype contrasts via t_test
    """
    connectivity_columns = [c for c in merged.columns if connectivity_col_contains in str(c)]
    if not connectivity_columns:
        raise ValueError(
            f"No connectivity columns found matching pattern '{connectivity_col_contains}'. "
            f"Tip: use --connectivity-col-contains to change the pattern."
        )

    pairwise_subtypes = ["AsymAD", "TypAD", "LowNFT"]
    pairwise_combinations = list(combinations(pairwise_subtypes, 2))

    # flatten all SNPs
    all_gwas_snps: Set[str] = {snp for snps in gwas_snp_groups.values() for snp in snps}

    formula_tmpl = "{conn} ~ {snp} * C(Subtype, Treatment('Control')) + age + sex + C(scanner) + C(scan_type)"

    results: List[dict] = []
    total = len(all_gwas_snps) * len(connectivity_columns)

    iterator = tqdm(all_gwas_snps, desc="SNPs", unit="snp") if progress else all_gwas_snps
    processed = 0

    for snp_col in iterator:
        for conn_col in connectivity_columns:
            processed += 1

            cols = [conn_col, snp_col, "Subtype", "age", "sex", "scanner", "scan_type"]
            # Drop NA like notebook
            try:
                subset = merged[cols].dropna()
            except KeyError:
                # If covariate columns missing, statsmodels can't run anyway
                continue

            if len(subset) < min_obs_for_model:
                continue
            if subset[snp_col].std() == 0 or subset[conn_col].std() == 0:
                continue
            if subset["Subtype"].nunique(dropna=True) < 2:
                continue

            row = {
                "SNP_Name": snp_col,
                "Connectivity_Name": conn_col,
                "R_squared": np.nan,
            }

            gwas_comp = [name for name, snps in gwas_snp_groups.items() if snp_col in snps]
            row["GWAS_Comparison"] = gwas_comp[0] if gwas_comp else "N/A"

            try:
                model = smf.ols(formula_tmpl.format(conn=conn_col, snp=snp_col), data=subset).fit()
                param_names = model.params.index.tolist()
                row["R_squared"] = float(model.rsquared)

                # subtype vs Control interaction terms
                for subtype in pairwise_subtypes:
                    param_name = f"{snp_col}:C(Subtype, Treatment('Control'))[T.{subtype}]"
                    if param_name in param_names:
                        row[f"P_{subtype}_vs_Control"] = float(model.pvalues[param_name])
                        row[f"Coeff_{subtype}_vs_Control"] = float(model.params[param_name])

                # custom contrasts: group1 - group2
                for g1, g2 in pairwise_combinations:
                    p1 = f"{snp_col}:C(Subtype, Treatment('Control'))[T.{g1}]"
                    p2 = f"{snp_col}:C(Subtype, Treatment('Control'))[T.{g2}]"
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
                # Mirror notebook behavior: silently skip failures
                pass

            if (not progress) and (processed % progress_every == 0):
                print(f"Processed {processed}/{total} combinations...")

    return pd.DataFrame(results)
