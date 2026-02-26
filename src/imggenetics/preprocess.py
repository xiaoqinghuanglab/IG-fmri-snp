from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from .utils import clean_column_name, normalize_subtype


def fisher_r_to_z(connectivity: pd.DataFrame) -> pd.DataFrame:
    """Apply Fisher r-to-z transform to float64 connectivity columns (excluding Subject_ID)."""
    cols = [c for c in connectivity.columns if c != "Subject_ID" and connectivity[c].dtype == "float64"]
    for col in cols:
        connectivity[col] = np.arctanh(np.clip(connectivity[col], -0.999999, 0.999999))
    return connectivity


def exclude_subject(
    genetic: pd.DataFrame,
    connectivity: pd.DataFrame,
    covariates: pd.DataFrame,
    subject_id: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not subject_id:
        return genetic, connectivity, covariates
    genetic = genetic[genetic["Subject_ID"] != subject_id]
    connectivity = connectivity[connectivity["Subject_ID"] != subject_id]
    covariates = covariates[covariates["Subject_ID"] != subject_id]
    return genetic, connectivity, covariates


def merge_for_analysis(genetic: pd.DataFrame, connectivity: pd.DataFrame, covariates: pd.DataFrame) -> pd.DataFrame:
    """Merge genetic + covariates + connectivity on Subject_ID (inner joins)."""
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
    """Prepare merged dataframe for statsmodels.

    Steps:
      - rename 'group' -> 'Subtype' if needed
      - normalize common subtype spellings (e.g., 'Typical AD' -> 'TypAD')
      - optionally drop subtypes not listed in subtype_categories
      - clean connectivity column names
      - cast Subtype, scanner, scan_type to categorical

    Returns:
      merged_prepped, map_original_to_cleaned_connectivity
    """
    if "group" in merged.columns and "Subtype" not in merged.columns:
        merged = merged.rename(columns={"group": "Subtype"})
    if "Subtype" not in merged.columns:
        raise KeyError("Missing 'Subtype' (or 'group') column in merged data.")

    merged["Subtype"] = merged["Subtype"].map(normalize_subtype)

    if drop_unknown_subtypes:
        merged = merged[merged["Subtype"].isin(subtype_categories)].copy()

    # Clean connectivity names (only those coming from the connectivity dataframe)
    original_connectivity_cols = [c for c in connectivity_df.columns if c != "Subject_ID"]
    name_map = {c: clean_column_name(c) for c in original_connectivity_cols}
    merged = merged.rename(columns=name_map)

    # Category casting
    if reference_subtype not in subtype_categories:
        raise ValueError(f"reference_subtype='{reference_subtype}' must be in subtype_categories={subtype_categories}")

    merged["Subtype"] = pd.Categorical(merged["Subtype"], categories=subtype_categories, ordered=False)
    if "scanner" in merged.columns:
        merged["scanner"] = pd.Categorical(merged["scanner"])
    if "scan_type" in merged.columns:
        merged["scan_type"] = pd.Categorical(merged["scan_type"])

    return merged, name_map
