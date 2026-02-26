#!/usr/bin/env python3
"""Export VEP rsID lists + per-GWAS result CSVs from imggenetics significant results.

This is a scriptified version of the post-processing done in the notebook (Cell-4),
but parameterized so paths can be supplied on the command line.

Inputs:
  --significant-csv : outputs/significant_result.csv produced by `imggenetics`
    must include columns: SNP_Name, GWAS_Comparison, Connectivity_Name, plus P/Coeff columns.

Outputs (under --out-dir):
  VEP_Input_Files/<GWAS_Comparison>_VEP.txt
  GWAS_Results_by_Comparison/<GWAS_Comparison>_results.csv

Example:
  python scripts/imggenetics/export_vep_inputs.py \
    --significant-csv outputs/significant_result.csv \
    --out-dir outputs/imggenetics_post
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


DEFAULT_GWAS_PAIRS = [
    "AsymAD_vs_Control",
    "Control_vs_TypAD",
    "AsymAD_vs_TypAD",
]


def parse_args():
    p = argparse.ArgumentParser(description="Export VEP input rsID lists and per-contrast CSVs.")
    p.add_argument("--significant-csv", type=Path, required=True, help="significant_result.csv from imggenetics.")
    p.add_argument("--out-dir", type=Path, required=True, help="Output folder.")
    p.add_argument("--gwas-pairs", nargs="+", default=DEFAULT_GWAS_PAIRS, help="Values in GWAS_Comparison column.")
    return p.parse_args()


def clean_rsid(snp_name: str) -> str:
    # remove trailing _A/_T/_C/_G
    return re.sub(r"_[ATCG]$", "", str(snp_name))


def main():
    args = parse_args()
    df = pd.read_csv(args.significant_csv)
    if "GWAS_Comparison" not in df.columns:
        raise SystemExit("Missing GWAS_Comparison column in significant CSV.")
    if "SNP_Name" not in df.columns:
        raise SystemExit("Missing SNP_Name column in significant CSV.")

    vep_dir = args.out_dir / "VEP_Input_Files"
    gwas_dir = args.out_dir / "GWAS_Results_by_Comparison"
    vep_dir.mkdir(parents=True, exist_ok=True)
    gwas_dir.mkdir(parents=True, exist_ok=True)

    for pair in args.gwas_pairs:
        sub = df[df["GWAS_Comparison"] == pair].copy()
        if sub.empty:
            print(f"[WARN] No rows for GWAS_Comparison='{pair}'. Skipping.")
            continue

        # VEP txt
        sub["SNP_Clean"] = sub["SNP_Name"].map(clean_rsid)
        rsids = sorted(sub["SNP_Clean"].dropna().astype(str).unique().tolist())
        vep_path = vep_dir / f"{pair}_VEP.txt"
        vep_path.write_text("\n".join(rsids) + "\n")

        # per-contrast CSV (keep all available columns)
        csv_path = gwas_dir / f"{pair}_results.csv"
        sub.to_csv(csv_path, index=False)

        print(f"[DONE] {pair}: rsIDs={len(rsids)} rows={len(sub)}")

    print(f"[DONE] Outputs written under: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
