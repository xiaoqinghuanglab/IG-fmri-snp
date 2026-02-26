#!/usr/bin/env python3
"""ADNI blood gene-expression preprocessing (probe collapse + covariate merge).

This script takes the raw ADNI expression export (wide matrix with header-like rows),
collapses multiple probes to one value per gene (keeping the probe with the highest
mean expression), and produces:

  1) a subject-by-gene expression matrix (CSV)
  2) a merged subject metadata table (CSV) that joins expression headers with your
     clinical covariates file

It is a parameterized version of the original pipeline so users can pass paths on the command line.

Example:
  python scripts/transcriptomics/adni_expression_preprocess.py \
    --expression-file ADNI_Gene_Expression_Profile.csv \
    --covariates-file covariates.csv \
    --out-matrix outputs/transcriptomics/adni_expression_matrix_cleaned.csv \
    --out-metadata outputs/transcriptomics/adni_metadata_merged.csv
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Preprocess ADNI gene expression file (probe collapse + covariate merge).")
    p.add_argument("--expression-file", type=Path, required=True, help="Raw ADNI expression CSV (as provided).")
    p.add_argument("--covariates-file", type=Path, required=True, help="Clinical covariates CSV for merging.")
    p.add_argument("--out-matrix", type=Path, required=True, help="Output subject-by-gene matrix CSV.")
    p.add_argument("--out-metadata", type=Path, required=True, help="Output merged metadata CSV.")
    p.add_argument("--subject-id-regex", type=str, default=r"(\d{3}_S_\d{4})", help="Regex to extract base subject ID.")
    return p.parse_args()


def main():
    args = parse_args()
    if not args.expression_file.exists():
        raise SystemExit(f"Expression file not found: {args.expression_file}")
    if not args.covariates_file.exists():
        raise SystemExit(f"Covariates file not found: {args.covariates_file}")

    print("[INFO] Loading expression file (may take time)...")
    data_full = pd.read_csv(args.expression_file, header=None, low_memory=False)

    # Metadata is typically in the top block (rows 0-7) with subject columns starting at col 3
    metadata_raw = data_full.iloc[0:8, 3:].T.reset_index(drop=True)
    metadata_raw.columns = [
        "Phase", "Visit", "SubjectID_Original", "purity_260_280",
        "purity_260_230", "rin", "affy_qc", "globin_qc",
    ]
    metadata_raw["SubjectID_Base"] = metadata_raw["SubjectID_Original"].astype(str).str.extract(args.subject_id_regex)

    # Expression block: row 8 onward
    gene_info = data_full.iloc[8:, :3].copy()
    gene_info.columns = ["ProbeSetID", "LocusLink", "GeneSymbol"]

    expression_values = data_full.iloc[8:, 3:].copy()
    expression_values = expression_values.apply(pd.to_numeric, errors="coerce")
    expression_values["GeneSymbol"] = gene_info["GeneSymbol"].values

    # Drop missing / placeholder symbols
    expression_values = expression_values.dropna(subset=["GeneSymbol"])
    expression_values = expression_values[expression_values["GeneSymbol"] != "---"]
    expression_values = expression_values[expression_values["GeneSymbol"].astype(str).str.strip() != ""]

    # Keep probe with highest mean per gene
    expression_values["Mean_Exp"] = expression_values.iloc[:, :-1].mean(axis=1)
    expression_values = expression_values.sort_values(by=["GeneSymbol", "Mean_Exp"])
    collapsed = expression_values.drop_duplicates(subset=["GeneSymbol"], keep="last")

    # Final matrix: index=GeneSymbol, columns=SubjectID_Original -> transpose to Subjects x Genes
    final_matrix = collapsed.set_index("GeneSymbol").drop(columns=["Mean_Exp"])
    final_matrix.columns = metadata_raw["SubjectID_Original"]
    final_matrix_T = final_matrix.T

    print("[INFO] Merging with covariates...")
    clinical_cov = pd.read_csv(args.covariates_file)

    # normalize common naming
    if "subject_id" in clinical_cov.columns and "SubjectID_Base" not in clinical_cov.columns:
        clinical_cov = clinical_cov.rename(columns={"subject_id": "SubjectID_Base"})

    merged_meta = pd.merge(metadata_raw, clinical_cov, on="SubjectID_Base", how="inner")
    common_subjects = merged_meta["SubjectID_Original"].values
    final_matrix_T = final_matrix_T.loc[common_subjects]

    args.out_matrix.parent.mkdir(parents=True, exist_ok=True)
    args.out_metadata.parent.mkdir(parents=True, exist_ok=True)
    final_matrix_T.to_csv(args.out_matrix)
    merged_meta.to_csv(args.out_metadata, index=False)

    max_val = float(pd.to_numeric(final_matrix_T.stack(), errors="coerce").max())
    print(f"[DONE] Matrix: {args.out_matrix} | shape={final_matrix_T.shape}")
    print(f"[DONE] Metadata: {args.out_metadata} | N={len(merged_meta)}")
    print(f"[CHECK] Max expression value: {max_val:.2f} (often <20 indicates already log-transformed)")


if __name__ == "__main__":
    main()
