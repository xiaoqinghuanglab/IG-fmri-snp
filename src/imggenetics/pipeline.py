from __future__ import annotations

from pathlib import Path
import pandas as pd

from .io import load_inputs, harmonize_subject_id
from .preprocess import fisher_r_to_z, exclude_subject, merge_for_analysis, prepare_for_model
from .model import build_gwas_snp_groups, run_ols_interaction_grid
from .fdr import apply_fdr_per_column, filter_significant_any_fdr


def run_pipeline(
    data_dir: Path,
    output_dir: Path,
    subject_to_exclude: str = "029_S_2395",
    connectivity_col_contains: str = "Connectivity",
    snp_prefix: str = "rs",
    min_obs_for_model: int = 10,
    progress: bool = True,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load
    genetic, connectivity, covariates, snp_metadata = load_inputs(data_dir)

    # Preprocess
    connectivity = fisher_r_to_z(connectivity)
    harmonize_subject_id(genetic, connectivity, covariates)
    genetic, connectivity, covariates = exclude_subject(genetic, connectivity, covariates, subject_to_exclude)

    merged = merge_for_analysis(genetic, connectivity, covariates)

    merged, _ = prepare_for_model(merged, connectivity_df=connectivity)

    # Build SNP groups and run regressions
    gwas_groups = build_gwas_snp_groups(merged, snp_metadata, snp_prefix=snp_prefix)

    results_df = run_ols_interaction_grid(
        merged,
        gwas_groups,
        connectivity_col_contains=connectivity_col_contains,
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
    sig_df = filter_significant_any_fdr(corrected_df)
    sig_path = output_dir / "significant_result.csv"
    sig_df.to_csv(sig_path, index=False)
