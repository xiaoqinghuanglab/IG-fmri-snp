#!/usr/bin/env python3
"""Export per-comparison VEP rsID lists from pairwise imggenetics results."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


RESULT_SUFFIX = "_significant_interaction_fdr05.csv"


def parse_args():
    parser = argparse.ArgumentParser(description="Export VEP rsID lists from pairwise SNP-connectivity results.")
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Output directory produced by imggenetics; comparison subfolders are scanned automatically.",
    )
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory for VEP files and copied CSVs.")
    return parser.parse_args()


def clean_rsid(value: str) -> str:
    return re.sub(r"_[ATCG]$", "", str(value))


def main():
    args = parse_args()
    results_dir = args.results_dir.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()

    result_files = sorted(results_dir.glob(f"*/*{RESULT_SUFFIX}"))
    if not result_files:
        raise SystemExit(f"No result files matching */*{RESULT_SUFFIX} were found under: {results_dir}")

    vep_dir = out_dir / "VEP_Input_Files"
    per_comparison_dir = out_dir / "GWAS_Results_by_Comparison"
    vep_dir.mkdir(parents=True, exist_ok=True)
    per_comparison_dir.mkdir(parents=True, exist_ok=True)

    for result_file in result_files:
        comparison = result_file.stem.removesuffix("_significant_interaction_fdr05")
        df = pd.read_csv(result_file)
        if "SNP" not in df.columns:
            raise SystemExit(f"Missing SNP column in: {result_file}")

        df = df.copy()
        df["SNP_Clean"] = df["SNP"].map(clean_rsid)
        rsids = sorted(df["SNP_Clean"].dropna().astype(str).unique().tolist())

        vep_path = vep_dir / f"{comparison}_VEP.txt"
        csv_path = per_comparison_dir / f"{comparison}_results.csv"
        vep_path.write_text("\n".join(rsids) + "\n")
        df.to_csv(csv_path, index=False)

        print(f"[DONE] {comparison}: rsIDs={len(rsids)} rows={len(df)}")

    print(f"[DONE] Outputs written under: {out_dir}")


if __name__ == "__main__":
    main()
