#!/usr/bin/env python3
"""Create (1) snp_metadata.csv and (2) genotype_matrix.csv from top-SNP lists.

This script is meant to sit at the boundary between PLINK and Python:
  - You provide a QC-passed PLINK prefix (--bfile)
  - You provide one or more *top SNP list* files (from extract_top_snps.py)
  - The script:
      (a) builds a combined metadata table (with case/control labels)
      (b) creates a union SNP list
      (c) calls PLINK to export numeric genotypes for those SNPs
      (d) formats a clean genotype_matrix.csv for downstream imaging-genetics

Example:
  python scripts/gwas/make_genotype_matrix.py \
    --plink plink \
    --bfile ADNI_Final_hg19_filtered \
    --top-list outputs/gwas/top1000_AsymAD_vs_Control.tsv AsymAD Control \
    --top-list outputs/gwas/top1000_AsymAD_vs_TypAD.tsv AsymAD TypAD \
    --top-list outputs/gwas/top1000_TypAD_vs_Control.tsv TypAD Control \
    --out-dir outputs/gwas/final
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd

from imggenetics.utils import normalize_subtype


def parse_args():
    p = argparse.ArgumentParser(description="Build snp_metadata.csv + genotype_matrix.csv from GWAS top lists.")
    p.add_argument("--plink", type=str, default="plink", help="PLINK executable (default: plink).")
    p.add_argument("--bfile", type=Path, required=True, help="QC-passed PLINK prefix (no extension).")
    p.add_argument(
        "--top-list",
        nargs=3,
        action="append",
        metavar=("TOP_TSV", "CASE", "CONTROL"),
        required=True,
        help="A top SNP TSV + labels. Repeat flag for multiple contrasts.",
    )
    p.add_argument("--out-dir", type=Path, required=True, help="Output directory.")
    p.add_argument("--top-n-col", type=str, default="SNP", help="Column in TOP_TSV containing rsID.")
    return p.parse_args()


def run(cmd: list[str]):
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    args = parse_args()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Build snp_metadata.csv
    rows = []
    for top_path, case, control in args.top_list:
        top_path = Path(top_path).expanduser().resolve()
        df = pd.read_csv(top_path, sep="\t")
        if args.top_n_col not in df.columns:
            raise SystemExit(f"Missing column '{args.top_n_col}' in {top_path}. Found: {list(df.columns)}")

        case_n = normalize_subtype(case)
        control_n = normalize_subtype(control)

        sub = df.copy()
        sub = sub.rename(columns={args.top_n_col: "SNP"})
        sub["Case"] = case_n
        sub["Control"] = control_n
        rows.append(sub[["SNP"] + [c for c in sub.columns if c not in {"SNP"}]])

    meta = pd.concat(rows, ignore_index=True)
    meta_out = out_dir / "snp_metadata.csv"
    meta.to_csv(meta_out, index=False)

    # 2) Union SNPs and export genotypes via PLINK
    union_snps = sorted(meta["SNP"].astype(str).unique().tolist())
    union_path = out_dir / "unique_snps_to_extract.txt"
    union_path.write_text("\n".join(union_snps) + "\n")

    plink_out_prefix = out_dir / "genotypes_union"
    run([
        args.plink,
        "--bfile", str(args.bfile),
        "--extract", str(union_path),
        "--recode", "A",
        "--out", str(plink_out_prefix),
    ])

    raw_path = Path(str(plink_out_prefix) + ".raw")
    if not raw_path.exists():
        raise SystemExit(f"PLINK .raw output not found: {raw_path}")

    raw = pd.read_csv(raw_path, delim_whitespace=True, dtype=str)
    # Standard columns: FID IID PAT MAT SEX PHENOTYPE + SNP columns
    if "IID" not in raw.columns:
        raise SystemExit(f"Unexpected .raw format. Columns: {list(raw.columns)[:20]}")

    geno = raw.copy()
    geno = geno.rename(columns={"IID": "Subject_ID"})
    keep_cols = ["Subject_ID"] + [c for c in geno.columns if c.startswith("rs")]
    geno = geno[keep_cols]

    # Convert genotype counts to numeric where possible
    for c in geno.columns:
        if c == "Subject_ID":
            continue
        geno[c] = pd.to_numeric(geno[c], errors="coerce")

    geno_out = out_dir / "genotype_matrix.csv"
    geno.to_csv(geno_out, index=False)

    print(f"[DONE] snp_metadata: {meta_out}")
    print(f"[DONE] genotype_matrix: {geno_out}")

if __name__ == "__main__":
    main()
