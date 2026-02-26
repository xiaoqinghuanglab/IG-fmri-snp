"""Multiple-testing correction utilities.

The pipeline generates many p-values (one per SNP × edge test, across several
subtype contrasts). We apply Benjamini–Hochberg (BH) FDR correction separately
for each p-value column.

Design choice:
  - Separate correction per column matches how the notebook reported contrast-
    specific significance (each contrast is its own family of tests).
"""

from __future__ import annotations

import pandas as pd
from statsmodels.stats import multitest


def apply_fdr_per_column(
    results_df: pd.DataFrame,
    alpha: float = 0.05,
    method: str = "fdr_bh",
) -> pd.DataFrame:
    """Apply FDR correction independently for each p-value column.

    Parameters
    ----------
    results_df:
        DataFrame produced by `run_ols_interaction_grid`.
    alpha:
        Target FDR level used by statsmodels (kept for completeness). We store the
        corrected p-values regardless of whether they pass alpha.
    method:
        statsmodels multipletests method; default 'fdr_bh'.

    Returns
    -------
    pd.DataFrame
        Copy of results_df with additional columns named:
            P_FDR<rest_of_original_name>
        Example:
            P_TypAD_vs_Control  ->  P_FDR_TypAD_vs_Control
    """

    corrected = results_df.copy()

    # By convention, the pipeline names p-value columns as:
    #   P_<group1>_vs_<group2>
    p_cols = [c for c in corrected.columns if str(c).startswith("P_") and "_vs_" in str(c)]

    for p_col in p_cols:
        valid = corrected[p_col].dropna()
        if len(valid) == 0:
            continue

        # multipletests returns (reject, pvals_corrected, ...)
        _, pvals_corr, _, _ = multitest.multipletests(pvals=valid.values, alpha=alpha, method=method)

        # Write back only where p-values exist (preserve NaNs).
        idx = valid.index
        corrected.loc[idx, f"P_FDR{str(p_col)[1:]}"] = pvals_corr

    return corrected


def filter_significant_any_fdr(corrected_df: pd.DataFrame, fdr_threshold: float = 0.05) -> pd.DataFrame:
    """Filter rows where *any* FDR-corrected p-value is below a threshold."""

    fdr_cols = [c for c in corrected_df.columns if str(c).startswith("P_FDR")]
    if not fdr_cols:
        # If no FDR columns exist, return an empty table (explicitly).
        return corrected_df.iloc[0:0].copy()

    mask = (corrected_df[fdr_cols] < fdr_threshold).any(axis=1)
    return corrected_df.loc[mask].copy()
