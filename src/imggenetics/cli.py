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
    p.add_argument("--data-dir", type=Path, default=Path("data"), help="Folder containing the 4 required CSV files.")
    p.add_argument("--output-dir", type=Path, default=Path("outputs"), help="Folder where outputs will be written.")
    p.add_argument("--subject-to-exclude", type=str, default="029_S_2395", help="Subject_ID to exclude from all datasets.")
    p.add_argument("--connectivity-col-contains", type=str, default="Connectivity",
                   help="Substring used to identify connectivity feature columns.")
    p.add_argument("--snp-prefix", type=str, default="rs", help="Prefix used to identify SNP columns (default: 'rs').")
    p.add_argument("--min-obs-for-model", type=int, default=10, help="Minimum N required to fit a model.")
    p.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bars.")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)

    run_pipeline(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        subject_to_exclude=args.subject_to_exclude,
        connectivity_col_contains=args.connectivity_col_contains,
        snp_prefix=args.snp_prefix,
        min_obs_for_model=args.min_obs_for_model,
        progress=not args.no_progress,
    )
    print(f"Done. Outputs written to: {args.output_dir.resolve()}")
    return 0
