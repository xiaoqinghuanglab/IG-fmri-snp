#!/usr/bin/env python3
"""Create PLINK phenotype + covariate files for subtype contrasts.

Input: a master CSV with one row per subject containing:
  - Subject_ID (or IID) and optionally FID
  - subtype label column (default: Subtype)
  - covariates: Age, Sex, APOE4, PC1..PC5

Output (per contrast):
  - pheno_<TAG>.txt  (FID IID PHENO where PHENO: 2=case, 1=control)
  - covar_<TAG>.txt  (FID IID Age Sex APOE4 PC1..PC5)

Example:
  python scripts/gwas/prepare_pheno_covar.py \
    --clinical clinical_data.csv \
    --out-dir outputs/gwas/pheno_covar \
    --subtype-col Subtypes \
    --id-col Subject_ID \
    --contrasts AsymAD:Control AsymAD:TypAD TypAD:Control
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from imggenetics.utils import normalize_subtype


def parse_args():
    p = argparse.ArgumentParser(description="Prepare PLINK pheno/covar files for subtype GWAS contrasts.")
    p.add_argument("--clinical", type=Path, required=True, help="Master clinical CSV.")
    p.add_argument("--out-dir", type=Path, required=True, help="Output folder.")
    p.add_argument("--subtype-col", type=str, default="Subtype", help="Column containing subtype labels.")
    p.add_argument("--id-col", type=str, default="Subject_ID", help="Subject identifier column (IID).")
    p.add_argument("--fid-col", type=str, default=None, help="Family ID column. If omitted, uses IID as FID.")
    p.add_argument("--age-col", type=str, default="Age", help="Age column.")
    p.add_argument("--sex-col", type=str, default="Sex", help="Sex column (numeric or string).")
    p.add_argument("--apoe4-col", type=str, default="APOE4", help="APOE4 carrier column.")
    p.add_argument("--pc-cols", nargs="+", default=["PC1", "PC2", "PC3", "PC4", "PC5"], help="PC covariate columns.")
    p.add_argument(
        "--contrasts",
        nargs="+",
        default=["AsymAD:Control", "AsymAD:TypAD", "TypAD:Control"],
        help="Contrasts as CASE:CONTROL (e.g., AsymAD:TypAD).",
    )
    return p.parse_args()


def encode_sex(s):
    # PLINK typically expects 1/2 or allows any numeric if --allow-no-sex
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce")
    x = s.astype(str).str.strip().str.lower()
    return x.map({"m": 1, "male": 1, "1": 1, "f": 2, "female": 2, "2": 2})


def main():
    args = parse_args()
    df = pd.read_csv(args.clinical)

    if args.id_col not in df.columns:
        raise SystemExit(f"Missing id column: {args.id_col}. Found: {list(df.columns)}")
    if args.subtype_col not in df.columns:
        raise SystemExit(f"Missing subtype column: {args.subtype_col}. Found: {list(df.columns)}")

    df = df.copy()
    df["IID"] = df[args.id_col].astype(str).str.strip()
    df["FID"] = df[args.fid_col].astype(str).str.strip() if args.fid_col and args.fid_col in df.columns else df["IID"]

    df["Subtype"] = df[args.subtype_col].map(normalize_subtype)
    df["Age"] = pd.to_numeric(df[args.age_col], errors="coerce")
    df["Sex"] = encode_sex(df[args.sex_col])
    df["APOE4"] = pd.to_numeric(df[args.apoe4_col], errors="coerce")

    for pc in args.pc_cols:
        if pc not in df.columns:
            raise SystemExit(f"Missing PC column: {pc}")
        df[pc] = pd.to_numeric(df[pc], errors="coerce")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    for c in args.contrasts:
        case, control = c.split(":")
        case = normalize_subtype(case)
        control = normalize_subtype(control)
        tag = f"{case}_vs_{control}"

        sub = df[df["Subtype"].isin([case, control])].copy()
        sub = sub.dropna(subset=["Age", "Sex", "APOE4"] + args.pc_cols)

        if sub.empty:
            print(f"[WARN] No rows for {tag}. Skipping.")
            continue

        sub["PHENO"] = (sub["Subtype"] == case).map({True: 2, False: 1}).astype(int)

        pheno = sub[["FID", "IID", "PHENO"]]
        covar = sub[["FID", "IID", "Age", "Sex", "APOE4"] + args.pc_cols]

        pheno_path = args.out_dir / f"pheno_{tag}.txt"
        covar_path = args.out_dir / f"covar_{tag}.txt"
        pheno.to_csv(pheno_path, sep="\t", index=False)
        covar.to_csv(covar_path, sep="\t", index=False)

        print(f"[DONE] {tag}: pheno={pheno_path.name} covar={covar_path.name} N={len(sub)}")

if __name__ == "__main__":
    main()
