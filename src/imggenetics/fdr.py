from __future__ import annotations

from typing import List

import pandas as pd
from statsmodels.stats import multitest


def apply_fdr_per_column(results_df: pd.DataFrame, alpha: float = 0.05, method: str = "fdr_bh") -> pd.DataFrame:
    """Apply BH-FDR correction separately for each p-value column.

    For each p-value column P_*_vs_*: writes P_FDR* column with corrected values.
    """
    corrected = results_df.copy()
    p_cols = [c for c in corrected.columns if str(c).startswith("P_") and "_vs_" in str(c)]

    for p_col in p_cols:
        valid = corrected[p_col].dropna()
        if len(valid) == 0:
            continue
        idx = valid.index
        _, pvals_corr, _, _ = multitest.multipletests(pvals=valid.values, alpha=alpha, method=method)
        corrected.loc[idx, f"P_FDR{str(p_col)[1:]}"] = pvals_corr

    return corrected


def filter_significant_any_fdr(corrected_df: pd.DataFrame, fdr_threshold: float = 0.05) -> pd.DataFrame:
    """Keep rows where ANY P_FDR* column is < threshold."""
    fdr_cols = [c for c in corrected_df.columns if str(c).startswith("P_FDR")]
    if not fdr_cols:
        return corrected_df.iloc[0:0].copy()
    mask = (corrected_df[fdr_cols] < fdr_threshold).any(axis=1)
    return corrected_df.loc[mask].copy()
