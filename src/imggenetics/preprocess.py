from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

from .utils import clean_column_name

# preprocess.py

import pandas as pd

def remove_low_nft(covariate_df):
    """
    Remove lowNFT group subjects.
    """
    print("Original subjects:", len(covariate_df))

    filtered_df = covariate_df[covariate_df["group"] != "lowNFT"].copy()

    print("After removing lowNFT:", len(filtered_df))
    return filtered_df


def generate_comparisons():
    """
    Define only the 3 valid comparisons.
    """
    comparisons = [
        ("AsymAD", "TypicalAD"),
        ("TypicalAD", "CN"),
        ("AsymAD", "CN")
    ]
    return comparisons


def fisher_r_to_z(connectivity: pd.DataFrame) -> pd.DataFrame:
    """
    Apply Fisher r-to-z transform to float64 connectivity columns (excluding Subject_ID).
    Mirrors the notebook logic.
    """
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
    """
    Merge genetic + covariates + connectivity on Subject_ID (inner joins), like the notebook.
    """
    merged = pd.merge(genetic, covariates, on="Subject_ID", how="inner")
    merged = pd.merge(merged, connectivity, on="Subject_ID", how="inner")
    return merged


def prepare_for_model(
    merged: pd.DataFrame,
    connectivity_df: pd.DataFrame,
    subtype_categories: List[str] = ["Control", "TypAD", "AsymAD", "LowNFT"],
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Prepare merged dataframe for statsmodels:
      - rename 'group' -> 'Subtype' if present
      - clean connectivity column names (but keep Subject_ID)
      - cast Subtype, scanner, scan_type to categorical

    Returns:
      merged_prepped, map_original_to_cleaned_connectivity
    """
    if "group" in merged.columns and "Subtype" not in merged.columns:
        merged = merged.rename(columns={"group": "Subtype"})

    if "Subtype" not in merged.columns:
        raise KeyError("Missing 'Subtype' (or 'group') column in merged data.")

    # Clean connectivity names (only those coming from the connectivity dataframe)
    original_connectivity_cols = [c for c in connectivity_df.columns if c != "Subject_ID"]
    name_map = {c: clean_column_name(c) for c in original_connectivity_cols}
    merged = merged.rename(columns=name_map)

    # Category casting
    merged["Subtype"] = pd.Categorical(merged["Subtype"], categories=subtype_categories, ordered=False)
    if "scanner" in merged.columns:
        merged["scanner"] = pd.Categorical(merged["scanner"])
    if "scan_type" in merged.columns:
        merged["scan_type"] = pd.Categorical(merged["scan_type"])

    return merged, name_map
