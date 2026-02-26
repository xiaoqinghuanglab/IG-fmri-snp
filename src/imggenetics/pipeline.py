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
    """End-to-end SNP×Subtype interaction pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load inputs
    if all(p is not None for p in [genetic_path, connectivity_path, covariates_path, snp_metadata_path]):
        genetic, connectivity, covariates, snp_metadata = load_inputs_from_paths(
            genetic_path=genetic_path, connectivity_path=connectivity_path,
            covariates_path=covariates_path, snp_metadata_path=snp_metadata_path,
        )
    else:
        genetic, connectivity, covariates, snp_metadata = load_inputs(data_dir)


    # Preprocess
    connectivity = fisher_r_to_z(connectivity)
    harmonize_subject_id(genetic, connectivity, covariates)
    genetic, connectivity, covariates = exclude_subject(genetic, connectivity, covariates, subject_to_exclude)

    merged = merge_for_analysis(genetic, connectivity, covariates)

    if subtype_categories is None:
        subtype_categories = ["Control", "TypAD", "AsymAD"]

    merged, _ = prepare_for_model(
        merged,
        connectivity_df=connectivity,
        subtype_categories=subtype_categories,
        reference_subtype=reference_subtype,
        drop_unknown_subtypes=drop_unknown_subtypes,
    )

    # Build SNP groups + run regressions
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

    # Save raw results
    raw_path = output_dir / "all_regression_result.csv"
    results_df.to_csv(raw_path, index=False)

    # FDR
    corrected_df = apply_fdr_per_column(results_df)
    corrected_path = output_dir / "all_regression_result_FDR.csv"
    corrected_df.to_csv(corrected_path, index=False)

    # Significant
    sig_df = filter_significant_any_fdr(corrected_df, fdr_threshold=fdr_threshold)
    sig_path = output_dir / "significant_result.csv"
    sig_df.to_csv(sig_path, index=False)
