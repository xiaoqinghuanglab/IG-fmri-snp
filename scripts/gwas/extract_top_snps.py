#!/usr/bin/env python3
"""Extract top-N SNPs from a PLINK *.assoc.logistic file (TEST=ADD).

Example:
  python scripts/gwas/extract_top_snps.py \
    --assoc gwas_AsymAD_vs_Control.assoc.logistic \
    --top-n 1000 \
    --out outputs/gwas/top1000_AsymAD_vs_Control.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Extract top-N SNPs from PLINK logistic output.")
    p.add_argument("--assoc", type=Path, required=True, help="PLINK .assoc.logistic file.")
    p.add_argument("--top-n", type=int, default=1000, help="Number of SNPs to keep.")
    p.add_argument("--out", type=Path, required=True, help="Output TSV.")
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.assoc, delim_whitespace=True, dtype=str)
    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"].copy()

    # Coerce numeric columns
    for c in ["P", "OR", "BP", "CHR"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["SNP", "P"]).sort_values("P", ascending=True)
    out = df.head(args.top_n).copy()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[DONE] Wrote {args.out} | rows={len(out)}")

if __name__ == "__main__":
    main()
