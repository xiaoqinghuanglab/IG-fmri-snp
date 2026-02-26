#!/usr/bin/env python3
"""ANMerge blood differential expression (covariate-adjusted OLS).

This is a parameterized version of the original ANMerge DEG script.
It runs OLS per gene with a group indicator + covariates, then BH-FDR.

Input CSV requirements:
  - subtype column (default: subtypes)
  - covariates (default: Age, sex)
  - remaining columns are treated as gene expression features

Example:
  python scripts/transcriptomics/anmerge_deg_covariates.py \
    --input ANMerge_Unified_Ready_for_DEA.csv \
    --out-dir outputs/transcriptomics/deg_anmerge
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from imggenetics.utils import normalize_subtype


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    m = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order]
    ranks = np.arange(1, m + 1)
    q = (ranked * m) / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    return out


def parse_args():
    p = argparse.ArgumentParser(description="ANMerge DEG with covariates (OLS per gene).")
    p.add_argument("--input", type=Path, required=True, help="ANMerge unified expression CSV.")
    p.add_argument("--out-dir", type=Path, required=True, help="Output folder.")
    p.add_argument("--dx-col", type=str, default="subtypes", help="Subtype column.")
    p.add_argument("--covariates", nargs="+", default=["Age", "sex"], help="Covariate column names.")
    p.add_argument(
        "--comparisons",
        nargs="+",
        default=["AsymAD:TypAD", "AsymAD:Control", "TypAD:Control"],
        help="Comparisons as TEST:REF.",
    )
    return p.parse_args()


def build_design(sub: pd.DataFrame, group1: str, dx_col: str, covariates: list[str]):
    group_ind = (sub[dx_col] == group1).astype(int)
    cov_df = sub[covariates].copy()
    cov_df = pd.get_dummies(cov_df, drop_first=True)
    X_df = pd.concat([group_ind.rename("__group__"), cov_df], axis=1)
    X = np.column_stack([np.ones(len(X_df), dtype=float), X_df.to_numpy(dtype=float)])
    row_mask = np.isfinite(X).all(axis=1)
    return X, row_mask


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    df[args.dx_col] = df[args.dx_col].map(normalize_subtype)

    clinical_cols = ["Subject_Id", args.dx_col] + args.covariates
    gene_cols = [c for c in df.columns if c not in clinical_cols]
    if not gene_cols:
        raise SystemExit("No gene columns detected. Check input file format.")

    for comp in args.comparisons:
        g1, g2 = comp.split(":")
        g1 = normalize_subtype(g1)
        g2 = normalize_subtype(g2)

        sub = df[df[args.dx_col].isin([g1, g2])].copy()
        if sub.empty:
            print(f"[WARN] no samples for {g1} vs {g2}")
            continue

        X, base_mask = build_design(sub, g1, args.dx_col, args.covariates)
        p = X.shape[1]
        results = []

        for gene in gene_cols:
            y = pd.to_numeric(sub[gene], errors="coerce").to_numpy(dtype=float)
            mask = base_mask & np.isfinite(y)
            n = int(mask.sum())
            if n <= p + 1:
                continue

            Xm = X[mask, :]
            ym = y[mask]
            if np.linalg.matrix_rank(Xm) < p:
                continue

            beta, _, _, _ = np.linalg.lstsq(Xm, ym, rcond=None)
            resid = ym - (Xm @ beta)
            df_resid = n - p
            if df_resid <= 0:
                continue

            s2 = (resid @ resid) / df_resid
            XtX_inv = np.linalg.inv(Xm.T @ Xm)
            se = np.sqrt(np.diag(s2 * XtX_inv))
            if se[1] == 0 or not np.isfinite(se[1]):
                continue

            t_stat = beta[1] / se[1]
            p_val = 2.0 * stats.t.sf(np.abs(t_stat), df=df_resid)

            results.append({
                "Gene": gene,
                "Beta_Group": float(beta[1]),
                "P_Value": float(p_val),
                "T_Stat": float(t_stat),
                "N_Used": int(n),
            })

        if not results:
            print(f"[WARN] no results for {g1} vs {g2}")
            continue

        res = pd.DataFrame(results).dropna(subset=["P_Value"])
        res["FDR"] = bh_fdr(res["P_Value"].to_numpy())
        res = res.sort_values("P_Value")

        out_path = args.out_dir / f"{g1}_vs_{g2}_CovAdj_OLS_DEG.csv"
        res.to_csv(out_path, index=False)
        print(f"[DONE] {out_path} | sig(FDR<0.05)={int((res['FDR']<0.05).sum())}")

    print(f"\n[DONE] ANMerge DEG finished: {args.out_dir}")

if __name__ == "__main__":
    main()
