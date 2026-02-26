"""High-level orchestration for the imaging-genetics OLS pipeline.

`run_pipeline()` stitches together:
  1) loading / subject-ID harmonization
  2) connectivity preprocessing (Fisher r-to-z)
  3) merging inputs into one modeling table
  4) building GWAS-derived SNP groups
  5) fitting SNP×Subtype interaction models for each SNP×edge
  6) multiple-testing correction (BH-FDR)
  7) exporting raw, FDR-corrected, and significant result tables

The goal is to keep the CLI thin and keep the scientific logic here.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd

from .fdr import apply_fdr_per_column, filter_significant_any_fdr
from .io import harmonize_subject_id, load_inputs, load_inputs_from_paths
from .model import build_gwas_snp_groups, run_ols_interaction_grid
from .preprocess import exclude_subject, fisher_r_to_z, merge_for_analysis, prepare_for_model


def run_pipeline(
    data_dir: Path,
    output_dir: Path,
    genetic_path: Optional[Path] = None,
    connectivity_path: Optional[Path] = None,
    covariates_path: Optional[Path] = None,
    snp_metadata_path: Optional[Path] = None,
    subject_to_exclude: str = "029_S_2395",
    connectivity_col_contains: str = "Connectivity",
    snp_prefix: str = "rs",
    min_obs_for_model: int = 10,
    progress: bool = True,
    subtype_categories: Optional[List[str]] = None,
    reference_subtype: str = "Control",
    drop_unknown_subtypes: bool = True,
    fdr_threshold: float = 0.05,
) -> None:
    """Run the full end-to-end pipeline.

    Parameters
    ----------
    data_dir:
        Folder containing default input filenames (see README).
    output_dir:
        Where to write result CSV files.
    genetic_path, connectivity_path, covariates_path, snp_metadata_path:
        Optional explicit input paths. If all four are provided, they override
        the data_dir lookup.
    subject_to_exclude:
        Optional subject ID to drop across all modalities.
    connectivity_col_contains:
        Pattern used to detect connectivity edge columns in the merged table.
    snp_prefix:
        Prefix to detect SNP columns.
    min_obs_for_model:
        Minimum complete-case N per SNP×edge test.
    progress:
        Whether to show a progress bar.
    subtype_categories:
        Subtype labels to include (and their order for categorical encoding).
    reference_subtype:
        Reference level used for treatment coding in statsmodels.
    drop_unknown_subtypes:
        If True, drop rows with subtype labels not present in subtype_categories.
    fdr_threshold:
        Threshold used when exporting significant_result.csv.
    """

    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # 1) Load inputs (either by explicit paths or default filenames)
    # ------------------------------------------------------------
    if all(p is not None for p in [genetic_path, connectivity_path, covariates_path, snp_metadata_path]):
        genetic, connectivity, covariates, snp_metadata = load_inputs_from_paths(
            genetic_path=genetic_path,
            connectivity_path=connectivity_path,
            covariates_path=covariates_path,
            snp_metadata_path=snp_metadata_path,
        )
    else:
        genetic, connectivity, covariates, snp_metadata = load_inputs(data_dir)

    # ------------------------------------------------------------
    # 2) Preprocess + harmonize join keys
    # ------------------------------------------------------------
    # Fisher r-to-z is applied to connectivity values to stabilize variance.
    connectivity = fisher_r_to_z(connectivity)

    # Standardize Subject_ID column header across all inputs.
    harmonize_subject_id(genetic, connectivity, covariates)

    # Optionally remove a known problematic subject.
    genetic, connectivity, covariates = exclude_subject(genetic, connectivity, covariates, subject_to_exclude)

    # ------------------------------------------------------------
    # 3) Merge inputs into one modeling table
    # ------------------------------------------------------------
    merged = merge_for_analysis(genetic, connectivity, covariates)

    # Default subtype categories if not provided.
    if subtype_categories is None:
        subtype_categories = ["Control", "TypAD", "AsymAD"]

    # Prepare categorical encoding and sanitize connectivity names.
    merged, _ = prepare_for_model(
        merged,
        connectivity_df=connectivity,
        subtype_categories=subtype_categories,
        reference_subtype=reference_subtype,
        drop_unknown_subtypes=drop_unknown_subtypes,
    )

    # ------------------------------------------------------------
    # 4) Build GWAS SNP groups + run regressions
    # ------------------------------------------------------------
    gwas_groups = build_gwas_snp_groups(merged, snp_metadata, snp_prefix=snp_prefix)

    results_df = run_ols_interaction_grid(
        merged,
        gwas_groups,
        connectivity_col_contains=connectivity_col_contains,
        reference_subtype=reference_subtype,
        subtype_categories=subtype_categories,
        min_obs_for_model=min_obs_for_model,
        progress=progress,
    )

    # ------------------------------------------------------------
    # 5) Save outputs
    # ------------------------------------------------------------
    # (a) Raw results
    raw_path = output_dir / "all_regression_result.csv"
    results_df.to_csv(raw_path, index=False)

    # (b) BH-FDR corrected p-values
    corrected_df = apply_fdr_per_column(results_df)
    corrected_path = output_dir / "all_regression_result_FDR.csv"
    corrected_df.to_csv(corrected_path, index=False)

    # (c) Significant subset (any contrast passes P_FDR < threshold)
    sig_df = filter_significant_any_fdr(corrected_df, fdr_threshold=fdr_threshold)
    sig_path = output_dir / "significant_result.csv"
    sig_df.to_csv(sig_path, index=False)
