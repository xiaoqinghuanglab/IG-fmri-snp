#!/usr/bin/env python3
"""ADNI blood differential expression (covariate-adjusted OLS).

Uses statsmodels OLS per gene:
  Expr ~ Group + Age + Sex + APOE4

Inputs:
  --matrix   CSV where rows=subjects and cols=genes (output of adni_expression_preprocess.py)
  --metadata CSV with at least: SubjectID_Original, Subtypes, Age, Gender(or Gender_Code), APOE4

Output:
  Writes one CSV per comparison to --out-dir, including BH-FDR.

Example:
  python scripts/transcriptomics/adni_deg_covariates.py \
    --matrix outputs/transcriptomics/adni_expression_matrix_cleaned.csv \
    --metadata outputs/transcriptomics/adni_metadata_merged.csv \
    --out-dir outputs/transcriptomics/deg_adni
"""

from __future__ import annotations

from typing import Optional

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from imggenetics.utils import normalize_subtype


def parse_args():
    p = argparse.ArgumentParser(description="ADNI DEG with covariates (OLS per gene).")
    p.add_argument("--matrix", type=Path, required=True, help="Subjects x genes CSV.")
    p.add_argument("--metadata", type=Path, required=True, help="Merged metadata CSV.")
    p.add_argument("--out-dir", type=Path, required=True, help="Output folder.")
    p.add_argument(
        "--comparisons",
        nargs="+",
        default=["AsymAD:TypAD", "AsymAD:Control", "TypAD:Control"],
        help="Comparisons as TEST:REF (e.g., AsymAD:Control).",
    )
    p.add_argument("--subtype-col", type=str, default="Subtypes", help="Subtype column name in metadata.")
    p.add_argument("--id-col", type=str, default="SubjectID_Original", help="Subject ID column in metadata.")
    p.add_argument("--age-col", type=str, default="Age", help="Age column name.")
    p.add_argument("--sex-col", type=str, default="Gender", help="Sex column name (M/F or numeric).")
    p.add_argument("--apoe4-col", type=str, default="APOE4", help="APOE4 column name.")
    return p.parse_args()


def encode_sex(series: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    return series.astype("category").cat.codes


def run_deg(test_group: str, ref_group: str, meta_df: pd.DataFrame, expr_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    print(f"\n[INFO] DEG: {test_group} vs {ref_group}")

    sub = meta_df[meta_df["Subtype_norm"].isin([test_group, ref_group])].copy()
    sub["Group_Code"] = np.where(sub["Subtype_norm"] == test_group, 1, 0)

    target_subs = sub["SubjectID"].values
    valid_subjects = [s for s in target_subs if s in expr_df.index]

    n_test = int((sub["Subtype_norm"] == test_group).isin([True]).sum())
    n_ref = int((sub["Subtype_norm"] == ref_group).isin([True]).sum())

    if len(valid_subjects) < 6:
        print("  [SKIP] Too few overlapping subjects.")
        return None

    current_meta = sub.set_index("SubjectID").loc[valid_subjects]
    current_expr = expr_df.loc[valid_subjects]

    X = current_meta[["Group_Code", "Age", "Sex_Code", "APOE4"]]
    X = sm.add_constant(X)
    Xv = X.values

    Y = current_expr.values
    genes = current_expr.columns

    results = []
    start = time.time()
    for i in range(len(genes)):
        y = Y[:, i]
        try:
            fit = sm.OLS(y, Xv).fit()
            results.append({
                "GeneSymbol": genes[i],
                "Beta_Group": float(fit.params[1]),
                "P_Value": float(fit.pvalues[1]),
            })
        except Exception:
            continue

        if i % 2000 == 0 and i > 0:
            sys.stdout.write(f"\r   progress {i}/{len(genes)} genes")
            sys.stdout.flush()

    print(f"\r   finished in {time.time()-start:.1f}s")

    res = pd.DataFrame(results).dropna(subset=["P_Value"]).sort_values("P_Value")
    if res.empty:
        return None

    res["FDR"] = multipletests(res["P_Value"], alpha=0.05, method="fdr_bh")[1]
    return res


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    expr = pd.read_csv(args.matrix, index_col=0)
    meta = pd.read_csv(args.metadata)

    # Normalize / rename IDs
    meta = meta.rename(columns={args.id_col: "SubjectID"})
    meta["Subtype_norm"] = meta[args.subtype_col].map(normalize_subtype)
    meta["Age"] = pd.to_numeric(meta[args.age_col], errors="coerce")
    meta["Sex_Code"] = encode_sex(meta[args.sex_col])
    meta["APOE4"] = pd.to_numeric(meta[args.apoe4_col], errors="coerce")

    meta = meta.dropna(subset=["SubjectID", "Subtype_norm", "Age", "Sex_Code", "APOE4"]).copy()

    for comp in args.comparisons:
        test, ref = comp.split(":")
        test = normalize_subtype(test)
        ref = normalize_subtype(ref)
        out = run_deg(test, ref, meta, expr)
        if out is None:
            continue
        out_path = args.out_dir / f"DEG_{test}_vs_{ref}_Covariates.csv"
        out.to_csv(out_path, index=False)
        print(f"[DONE] {out_path} | sig(FDR<0.05)={int((out['FDR']<0.05).sum())}")

    print("\n[DONE] All ADNI DEG comparisons finished.")

if __name__ == "__main__":
    main()
