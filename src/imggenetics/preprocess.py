"""Preprocessing helpers for imaging-genetics modeling.

These functions implement the lightweight transformations needed to:
  - put connectivity values on a scale suitable for linear modeling
  - merge genetic, imaging, and covariate tables
  - enforce consistent subtype labels and categorical encodings

They are intentionally data-agnostic: most path / file assumptions live in
`imggenetics.io` and the CLI.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from .utils import clean_column_name, normalize_subtype


def fisher_r_to_z(connectivity: pd.DataFrame) -> pd.DataFrame:
    """Apply Fisher r-to-z transform to connectivity columns.

    Why:
      - Many connectivity matrices start as correlations (bounded in [-1, 1]).
      - Fisher transform stabilizes variance and makes correlations more
        approximately normal for regression.

    Notes:
      - We clip correlations very slightly to avoid infinities in arctanh.
      - Only numeric columns are transformed; `Subject_ID` is excluded.
    """

    # Transform only numeric feature columns.
    cols = [c for c in connectivity.columns if c != "Subject_ID" and pd.api.types.is_numeric_dtype(connectivity[c])]
    for col in cols:
        connectivity[col] = np.arctanh(np.clip(connectivity[col].astype(float), -0.999999, 0.999999))
    return connectivity


def exclude_subject(
    genetic: pd.DataFrame,
    connectivity: pd.DataFrame,
    covariates: pd.DataFrame,
    subject_id: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Drop a single subject across all inputs.

    This mirrors the original notebook behavior where a problematic subject was
    removed prior to modeling.
    """

    if not subject_id:
        return genetic, connectivity, covariates

    genetic = genetic[genetic["Subject_ID"] != subject_id]
    connectivity = connectivity[connectivity["Subject_ID"] != subject_id]
    covariates = covariates[covariates["Subject_ID"] != subject_id]
    return genetic, connectivity, covariates


def merge_for_analysis(genetic: pd.DataFrame, connectivity: pd.DataFrame, covariates: pd.DataFrame) -> pd.DataFrame:
    """Inner-join genetic + covariates + connectivity on Subject_ID."""

    merged = pd.merge(genetic, covariates, on="Subject_ID", how="inner")
    merged = pd.merge(merged, connectivity, on="Subject_ID", how="inner")
    return merged


def prepare_for_model(
    merged: pd.DataFrame,
    connectivity_df: pd.DataFrame,
    subtype_categories: List[str],
    reference_subtype: str = "Control",
    drop_unknown_subtypes: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """Prepare a merged table for statsmodels OLS.

    Steps
    -----
    1) Ensure there is a Subtype column (accepts `group` as an alias).
    2) Normalize subtype spellings (e.g., 'Typical AD' -> 'TypAD').
    3) Optionally drop rows with subtypes not listed in subtype_categories.
    4) Sanitize connectivity column names to be formula-safe.
    5) Convert Subtype/scanner/scan_type to categorical dtypes.

    Returns
    -------
    merged_prepped:
        The cleaned table used for modeling.
    name_map:
        Mapping from original connectivity column names to cleaned names.
        This is useful if you want to map results back to the original labels.
    """

    # Some upstream scripts use `group` instead of `Subtype`.
    if "group" in merged.columns and "Subtype" not in merged.columns:
        merged = merged.rename(columns={"group": "Subtype"})

    if "Subtype" not in merged.columns:
        raise KeyError("Missing 'Subtype' (or 'group') column in merged data.")

    # Normalize subtype strings (handles spacing/casing differences).
    merged["Subtype"] = merged["Subtype"].map(normalize_subtype)

    if drop_unknown_subtypes:
        merged = merged[merged["Subtype"].isin(subtype_categories)].copy()

    # Clean connectivity names (only those coming from the connectivity dataframe).
    original_connectivity_cols = [c for c in connectivity_df.columns if c != "Subject_ID"]
    name_map = {c: clean_column_name(c) for c in original_connectivity_cols}
    merged = merged.rename(columns=name_map)

    # Category casting ensures statsmodels uses consistent treatment coding.
    if reference_subtype not in subtype_categories:
        raise ValueError(
            f"reference_subtype='{reference_subtype}' must be in subtype_categories={subtype_categories}"
        )

    merged["Subtype"] = pd.Categorical(merged["Subtype"], categories=subtype_categories, ordered=False)

    # These columns are optional; include if present.
    if "scanner" in merged.columns:
        merged["scanner"] = pd.Categorical(merged["scanner"])
    if "scan_type" in merged.columns:
        merged["scan_type"] = pd.Categorical(merged["scan_type"])

    return merged, name_map
