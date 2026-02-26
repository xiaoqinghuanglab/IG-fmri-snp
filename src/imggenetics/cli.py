from __future__ import annotations

import argparse
from pathlib import Path

from .pipeline import run_pipeline
from .utils import setup_logging


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="imggenetics",
        description="Automated SNP×Subtype OLS interaction pipeline for connectivity phenotypes.",
    )
    p.add_argument("--data-dir", type=Path, default=Path("data"), help="Folder containing the 4 required CSV files (default filenames).")

    p.add_argument("--genetic-file", type=Path, default=None, help="Optional path to genotype_matrix.csv (overrides --data-dir lookup).")
    p.add_argument("--connectivity-file", type=Path, default=None, help="Optional path to connectivity CSV (overrides --data-dir lookup).")
    p.add_argument("--covariates-file", type=Path, default=None, help="Optional path to subject_covariates.csv (overrides --data-dir lookup).")
    p.add_argument("--snp-metadata-file", type=Path, default=None, help="Optional path to snp_metadata.csv (overrides --data-dir lookup).")
    p.add_argument("--output-dir", type=Path, default=Path("outputs"), help="Folder where outputs will be written.")
    p.add_argument("--subject-to-exclude", type=str, default="029_S_2395", help="Subject_ID to exclude.")
    p.add_argument("--connectivity-col-contains", type=str, default="Connectivity", help="Substring used to identify connectivity feature columns.")
    p.add_argument("--snp-prefix", type=str, default="rs", help="Prefix used to identify SNP columns (default: 'rs').")
    p.add_argument("--min-obs-for-model", type=int, default=10, help="Minimum N required to fit a model.")
    p.add_argument("--fdr-threshold", type=float, default=0.05, help="Significance threshold applied to any P_FDR* column.")

    p.add_argument(
        "--subtypes",
        nargs="+",
        default=["Control", "TypAD", "AsymAD"],
        help="Subtype categories to include (default: Control TypAD AsymAD).",
    )
    p.add_argument(
        "--reference-subtype",
        type=str,
        default="Control",
        help="Reference subtype for Treatment coding (default: Control).",
    )
    p.add_argument(
        "--keep-unknown-subtypes",
        action="store_true",
        help="Do not drop rows with subtype labels not listed in --subtypes.",
    )

    p.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bars.")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)

    run_pipeline(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        genetic_path=args.genetic_file,
        connectivity_path=args.connectivity_file,
        covariates_path=args.covariates_file,
        snp_metadata_path=args.snp_metadata_file,
        subject_to_exclude=args.subject_to_exclude,
        connectivity_col_contains=args.connectivity_col_contains,
        snp_prefix=args.snp_prefix,
        min_obs_for_model=args.min_obs_for_model,
        progress=not args.no_progress,
        subtype_categories=args.subtypes,
        reference_subtype=args.reference_subtype,
        drop_unknown_subtypes=not args.keep_unknown_subtypes,
        fdr_threshold=args.fdr_threshold,
    )

    print(f"Done. Outputs written to: {args.output_dir.resolve()}")
    return 0
