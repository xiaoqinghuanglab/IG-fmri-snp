#!/usr/bin/env python3
"""
ANMerge B2 Prep (GeneExpression -> MRI) with Contrast-specific Gene Lists
========================================================================

Inputs in your folder:
- ANMerge_MRI_FS6.0_under_90.csv
- ANMerge_NewModel_NewLabel_Clinical_Info.csv
- Batch Removed ANMerge_blood_rna_gene_expr_removedbatch_02022024.csv
- *UniqueGenes.csv  (e.g., "AsymAD vs TypAD_UniqueGenes.csv")

This script will:
1) Filter clinical to ONLY: Asym AD, Typical AD, Control
2) Select MRI baseline per subject
3) Compute overlap:
   - Clinical ∩ Expression
   - Clinical ∩ Expression ∩ MRI baseline
4) Save prepared subject-level tables:
   - prepared_gene_expression_sample.csv   (clinical + optional MRI meta; no full expr)
   - prepared_expr_mri_merged_sample.csv   (clinical + MRI baseline; no full expr)
5) Save demographics tables (like your screenshot)
6) NEW: For each *UniqueGenes.csv:
   - load gene list from column "Unique Gene Name" (fallback: first column)
   - extract expression for merged cohort subjects (N=200)
   - write:
       selected_gene_expression_expr_mri_<TAG>_subject_by_gene.csv
       prepared_expr_mri_merged_with_selected_genes_<TAG>.csv
       gene_list_<TAG>_found.txt
       gene_list_<TAG>_missing.txt

Run:
  cd ~/Desktop/zAnMerge
  python anmerge_b2_gene_mri.py --data-dir .

Optional:
  python anmerge_b2_gene_mri.py --data-dir . --gene-files "AsymAD vs TypAD_UniqueGenes.csv" "AsymAD vs Control_UniqueGenes.csv"
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# scipy for demographics p-values
try:
    from scipy import stats
except ImportError:
    print("[ERROR] scipy is required. Install: pip install scipy")
    sys.exit(1)


# -----------------------------
# Config
# -----------------------------
CANON_SUBTYPES = ["Asym AD", "Typical AD", "Control"]
SUBTYPE_ORDER  = ["Asym AD", "Typical AD", "Control"]


# -----------------------------
# Utils
# -----------------------------
def _clean_id(x) -> str:
    return str(x).strip()

def _norm_colname(x: str) -> str:
    x = str(x).strip().lower()
    x = re.sub(r"[^a-z0-9]+", "", x)
    return x

def _to_numeric_safe(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.replace([np.inf, -np.inf, "inf", "-inf"], np.nan), errors="coerce")

def fmt_mean_sd(x: pd.Series) -> str:
    x = pd.to_numeric(x, errors="coerce")
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return ""
    return f"{x.mean():.2f} ± {x.std(ddof=1):.2f}"

def fmt_count_pct(n: int, N: int) -> str:
    if N == 0:
        return f"{n} (0.00%)"
    return f"{n} ({100.0 * n / N:.2f}%)"

def fmt_pval(p: float) -> str:
    if p is None or (isinstance(p, float) and (np.isnan(p) or not np.isfinite(p))):
        return ""
    return f"{p:.2E}"

def normalize_subtype(v: str) -> str:
    if v is None:
        return ""
    s = str(v).strip().lower()
    s2 = re.sub(r"[\s_\-]+", "", s)

    if s2 in {"control", "cn", "ctl", "cognitivenormal", "normal"}:
        return "Control"
    if s2 in {"asymad", "asymptomaticad"} or s.startswith("asym"):
        return "Asym AD"
    if s2 in {"typicalad"} or s.startswith("typical"):
        return "Typical AD"
    return ""

def detect_subject_id_column(df: pd.DataFrame) -> str:
    """Detect likely subject-ID column by values like DCR00001, KPO..., etc."""
    id_like = re.compile(r"^[A-Za-z]{2,}\d{3,}.*$")
    for c in df.columns:
        s = df[c].astype(str).str.strip()
        sample = s[s.notna()].head(200)
        if len(sample) == 0:
            continue
        frac = (sample.str.match(id_like)).mean()
        if frac > 0.80:
            return c
    return ""

def slugify_filename(name: str) -> str:
    """Make a safe tag from filename."""
    base = Path(name).stem
    base = base.replace("UniqueGenes", "")
    base = base.replace("__", "_")
    base = re.sub(r"\s+", "_", base.strip())
    base = re.sub(r"[^A-Za-z0-9_]+", "", base)
    base = re.sub(r"_+", "_", base).strip("_")
    return base or "geneset"


# -----------------------------
# Loaders
# -----------------------------
def get_expr_header(expr_path: Path):
    header = pd.read_csv(expr_path, nrows=0)
    cols = [str(c).strip() for c in header.columns]
    if not cols:
        raise ValueError("Expression file has no columns.")
    gene_col = cols[0]
    subj_cols = cols[1:]
    return gene_col, subj_cols

def get_expr_subject_ids(expr_path: Path) -> pd.Index:
    _, subj_cols = get_expr_header(expr_path)
    return pd.Index([_clean_id(c) for c in subj_cols])

def load_mri_baseline(mri_path: Path) -> pd.DataFrame:
    mri = pd.read_csv(mri_path)
    mri.columns = [str(c).strip() for c in mri.columns]

    id_col = detect_subject_id_column(mri)
    if not id_col:
        id_col = mri.columns[0]

    mri[id_col] = mri[id_col].map(_clean_id)
    mri = mri.set_index(id_col)
    mri.index.name = "SUBJ"

    if "Month" in mri.columns:
        mri["Month"] = _to_numeric_safe(mri["Month"])
        mri_reset = mri.reset_index().sort_values(["SUBJ", "Month"])
    elif "Visit" in mri.columns:
        mri["Visit"] = _to_numeric_safe(mri["Visit"])
        mri_reset = mri.reset_index().sort_values(["SUBJ", "Visit"])
    else:
        raise ValueError("MRI file must contain 'Month' or 'Visit'.")

    mri_base = (
        mri_reset.groupby("SUBJ", as_index=False)
        .first()
        .set_index("SUBJ")
    )

    # Clean common metadata if present
    for c in ["Age", "MMSE"]:
        if c in mri_base.columns:
            mri_base[c] = _to_numeric_safe(mri_base[c])
    for c in ["Sex", "Site", "Diagnosis", "ScanType", "APOE"]:
        if c in mri_base.columns:
            mri_base[c] = mri_base[c].astype(str).str.strip()

    return mri_base

def load_clinical(clin_path: Path) -> pd.DataFrame:
    clin_raw = pd.read_csv(clin_path)
    clin_raw.columns = [str(c).strip() for c in clin_raw.columns]

    id_col = detect_subject_id_column(clin_raw)
    if id_col:
        clin_raw[id_col] = clin_raw[id_col].map(_clean_id)
        clin = clin_raw.set_index(id_col)
    else:
        clin = pd.read_csv(clin_path, index_col=0)
        clin.columns = [str(c).strip() for c in clin.columns]

    clin.index = clin.index.map(_clean_id)
    clin.index.name = "SUBJ"

    norm_map = {_norm_colname(c): c for c in clin.columns}

    if "gender" in norm_map and "Sex" not in clin.columns:
        clin = clin.rename(columns={norm_map["gender"]: "Sex"})
    if "sex" in norm_map and "Sex" not in clin.columns:
        clin = clin.rename(columns={norm_map["sex"]: "Sex"})

    subtype_col = None
    for key in ["subtypes", "subtype", "newlabel", "label"]:
        if key in norm_map:
            subtype_col = norm_map[key]
            break
    if subtype_col is None:
        print("\n[DEBUG] Clinical columns found:\n", list(clin.columns))
        raise ValueError("Clinical subtype column not found (expected 'subtypes' or similar).")

    clin = clin.rename(columns={subtype_col: "Subtype"})

    if "Sex" in clin.columns:
        clin["Sex"] = clin["Sex"].astype(str).str.strip()
    if "Age" in clin.columns:
        clin["Age"] = _to_numeric_safe(clin["Age"])

    clin["Subtype_raw"] = clin["Subtype"].astype(str).str.strip()
    clin["Subtype"] = clin["Subtype_raw"].map(normalize_subtype)

    clin = clin[clin["Subtype"].isin(CANON_SUBTYPES)].copy()
    clin["Subtype"] = pd.Categorical(clin["Subtype"], categories=SUBTYPE_ORDER, ordered=True)

    return clin


# -----------------------------
# Gene list loader
# -----------------------------
def load_gene_list_csv(path: Path) -> list:
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]

    # Prefer "Unique Gene Name" (exact or normalized match)
    target = None
    for c in df.columns:
        if c == "Unique Gene Name":
            target = c
            break
    if target is None:
        # normalized search
        nmap = {_norm_colname(c): c for c in df.columns}
        if "uniquegenename" in nmap:
            target = nmap["uniquegenename"]

    # fallback: first column
    if target is None:
        target = df.columns[0]

    genes = (
        df[target]
        .astype(str)
        .str.strip()
        .replace({"nan": np.nan, "None": np.nan, "": np.nan})
        .dropna()
        .tolist()
    )
    # keep unique in order
    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


# -----------------------------
# Expression subset (memory-safe)
# -----------------------------
def load_expression_subset_subject_by_gene(expr_path: Path, subject_ids, genes, chunksize: int = 2000) -> pd.DataFrame:
    """
    Load (subjects x genes) for selected genes only, restricted to given subject IDs.
    Returns DataFrame index=SUBJ, columns=EXPR_<GENE>.
    """
    gene_col, subj_cols = get_expr_header(expr_path)

    subj_set = set(map(str, subject_ids))
    keep_subjects = [c for c in subj_cols if _clean_id(c) in subj_set]
    usecols = [gene_col] + keep_subjects

    genes = [g.strip() for g in genes if str(g).strip()]
    gene_set = set(genes)
    if not gene_set:
        raise ValueError("Empty gene set.")

    kept_rows = []
    for chunk in pd.read_csv(expr_path, usecols=usecols, chunksize=chunksize):
        chunk[gene_col] = chunk[gene_col].astype(str).str.strip()
        sub = chunk[chunk[gene_col].isin(gene_set)]
        if not sub.empty:
            kept_rows.append(sub)

    if not kept_rows:
        raise ValueError("None of the requested genes were found in the expression file.")

    expr_sub = pd.concat(kept_rows, axis=0, ignore_index=True).set_index(gene_col)
    expr_sub.columns = [_clean_id(c) for c in expr_sub.columns]

    expr_t = expr_sub.T  # subjects x genes
    expr_t = expr_t.apply(pd.to_numeric, errors="coerce")
    expr_t = expr_t.rename(columns={g: f"EXPR_{g}" for g in expr_t.columns})
    expr_t.index.name = "SUBJ"
    return expr_t


# -----------------------------
# Demographics table builder
# -----------------------------
def demographics_table(df: pd.DataFrame, group_col: str = "Subtype") -> pd.DataFrame:
    groups = SUBTYPE_ORDER
    cols = groups + ["Total", "p_value"]

    if df is None or len(df) == 0:
        out = pd.DataFrame([{
            "Metric": "N",
            **{g: "N=0" for g in groups},
            "Total": "N=0",
            "p_value": ""
        }])
        return out[["Metric"] + cols]

    df = df[df[group_col].isin(CANON_SUBTYPES)].copy()
    df[group_col] = pd.Categorical(df[group_col], categories=SUBTYPE_ORDER, ordered=True)

    out_rows = []

    # N row
    n_by = {g: int((df[group_col] == g).sum()) for g in groups}
    n_tot = int(df[group_col].notna().sum())
    out_rows.append({"Metric": "N", **{g: f"N={n_by[g]}" for g in groups}, "Total": f"N={n_tot}", "p_value": ""})

    # Age
    if "Age" in df.columns:
        row = {"Metric": "Age (Mean (SD))"}
        for g in groups:
            row[g] = fmt_mean_sd(df.loc[df[group_col] == g, "Age"])
        row["Total"] = fmt_mean_sd(df["Age"])

        arrays = []
        ok = True
        for g in groups:
            x = pd.to_numeric(df.loc[df[group_col] == g, "Age"], errors="coerce")
            x = x[np.isfinite(x)]
            if len(x) < 2:
                ok = False
            arrays.append(x)
        p = stats.f_oneway(*arrays).pvalue if ok else np.nan
        row["p_value"] = fmt_pval(p)
        out_rows.append(row)

    # Categorical vars (include if present)
    cat_vars = []
    if "Sex" in df.columns:
        cat_vars.append(("Sex (n%)", "Sex"))
    if "Site" in df.columns:
        cat_vars.append(("Site (n%)", "Site"))
    if "Diagnosis" in df.columns:
        cat_vars.append(("Diagnosis (n%)", "Diagnosis"))
    if "ScanType" in df.columns:
        cat_vars.append(("Scan Type (n%)", "ScanType"))
    if "APOE" in df.columns:
        cat_vars.append(("APOE (n%)", "APOE"))

    for label, var in cat_vars:
        df_tmp = df.copy()
        df_tmp[var] = df_tmp[var].astype(str).str.strip().replace({"nan": np.nan, "None": np.nan})

        ct = pd.crosstab(df_tmp[group_col], df_tmp[var])
        try:
            if ct.shape[0] == len(groups) and ct.shape[1] >= 2:
                p_cat = stats.chi2_contingency(ct.loc[groups].values)[1]
            else:
                p_cat = np.nan
        except Exception:
            p_cat = np.nan

        levels = df_tmp[var].value_counts(dropna=True).index.tolist()
        first = True
        for lev in levels:
            row = {"Metric": f"{label}: {lev}"}
            for g in groups:
                sub = df_tmp.loc[df_tmp[group_col] == g, var]
                N = int(sub.notna().sum())
                n = int((sub == lev).sum())
                row[g] = fmt_count_pct(n, N)

            sub_all = df_tmp[var]
            N_all = int(sub_all.notna().sum())
            n_all = int((sub_all == lev).sum())
            row["Total"] = fmt_count_pct(n_all, N_all)

            row["p_value"] = fmt_pval(p_cat) if first else ""
            first = False
            out_rows.append(row)

    out = pd.DataFrame(out_rows)
    return out[["Metric"] + cols]


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=str, default=".", help="Folder containing CSVs")
    ap.add_argument("--mri",  type=str, default="ANMerge_MRI_FS6.0_under_90.csv")
    ap.add_argument("--clin", type=str, default="ANMerge_NewModel_NewLabel_Clinical_Info.csv")
    ap.add_argument("--expr", type=str, default="Batch Removed ANMerge_blood_rna_gene_expr_removedbatch_02022024.csv")

    # NEW: gene list files (optional). If not provided, script auto-detects *UniqueGenes.csv
    ap.add_argument("--gene-files", nargs="*", default=None,
                    help="One or more *UniqueGenes.csv files. If omitted, auto-detects in data-dir.")

    # Performance
    ap.add_argument("--chunksize", type=int, default=2000, help="CSV chunksize for expression reading")
    args = ap.parse_args()

    data_dir = Path(args.data_dir).expanduser().resolve()
    mri_path  = data_dir / args.mri
    clin_path = data_dir / args.clin
    expr_path = data_dir / args.expr

    for p in [mri_path, clin_path, expr_path]:
        if not p.exists():
            print(f"[ERROR] Missing file: {p}")
            sys.exit(1)

    out_dir = data_dir / "results_b2"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load main tables
    print("[INFO] Loading clinical...")
    clin = load_clinical(clin_path)
    print("  Clinical N:", len(clin))
    print("  Subtype counts:\n", clin["Subtype"].value_counts(dropna=False))

    print("[INFO] Loading MRI baseline...")
    mri_base = load_mri_baseline(mri_path)
    print("  MRI baseline N:", len(mri_base))

    print("[INFO] Reading expression subject IDs (header only)...")
    expr_ids = get_expr_subject_ids(expr_path)
    print("  Expression subject IDs N:", len(expr_ids))

    # Cohort A: clinical ∩ expression
    ids_expr = clin.index.intersection(expr_ids)
    df_expr = clin.loc[ids_expr].copy()

    # attach some MRI metadata to gene-expression sample (left join)
    meta_cols = [c for c in ["Site", "Diagnosis", "ScanType", "APOE", "MMSE"] if c in mri_base.columns]
    if meta_cols:
        df_expr = df_expr.join(mri_base.loc[mri_base.index.intersection(ids_expr), meta_cols], how="left")

    print("[INFO] Gene Expression sample (Clinical ∩ Expr):", len(df_expr))

    # Cohort B: clinical ∩ expression ∩ MRI baseline
    ids_expr_mri = clin.index.intersection(expr_ids).intersection(mri_base.index)
    df_expr_mri = clin.loc[ids_expr_mri].join(mri_base.loc[ids_expr_mri], how="inner", rsuffix="_MRI")
    print("[INFO] Expr+MRI merged sample (Clinical ∩ Expr ∩ MRI):", len(df_expr_mri))

    # Save prepared subject-level tables (no full expr)
    prep_expr_path = out_dir / "prepared_gene_expression_sample.csv"
    df_expr.reset_index().to_csv(prep_expr_path, index=False)
    print("  Saved:", prep_expr_path)

    prep_expr_mri_path = out_dir / "prepared_expr_mri_merged_sample.csv"
    df_expr_mri.reset_index().to_csv(prep_expr_mri_path, index=False)
    print("  Saved:", prep_expr_mri_path)

    # Save ID lists
    (out_dir / "ids_gene_expression_sample.txt").write_text("\n".join(map(str, ids_expr)) + "\n")
    (out_dir / "ids_expr_mri_merged_sample.txt").write_text("\n".join(map(str, ids_expr_mri)) + "\n")

    # Demographics
    print("[INFO] Building demographics table: Gene Expression sample...")
    demo_expr = demographics_table(df_expr, group_col="Subtype")
    demo_expr.to_csv(out_dir / "demographics_gene_expression.csv", index=False)

    print("[INFO] Building demographics table: Expr+MRI merged sample...")
    demo_expr_mri = demographics_table(df_expr_mri, group_col="Subtype")
    demo_expr_mri.to_csv(out_dir / "demographics_expr_mri_merged.csv", index=False)

    # Optional Excel with both tables
    xlsx_path = out_dir / "demographics_tables.xlsx"
    try:
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            demo_expr.to_excel(writer, sheet_name="GeneExpression", index=False)
            demo_expr_mri.to_excel(writer, sheet_name="Expr_MRI_Merged", index=False)
        print("  Saved:", xlsx_path)
    except Exception as e:
        print("[WARN] Excel write failed (openpyxl missing?). CSVs are saved. Error:", str(e))

    # -----------------------------
    # NEW: Process gene list CSVs
    # -----------------------------
    if args.gene_files is None:
        gene_files = sorted(data_dir.glob("*UniqueGenes.csv"))
    else:
        gene_files = [data_dir / f for f in args.gene_files]

    if not gene_files:
        print("[INFO] No *UniqueGenes.csv provided/found. Done.")
        return

    # Load all gene lists
    genesets = {}
    for gf in gene_files:
        if not gf.exists():
            print(f"[WARN] Gene file not found (skipping): {gf}")
            continue
        tag = slugify_filename(gf.name)
        genes = load_gene_list_csv(gf)
        genesets[tag] = genes
        print(f"[INFO] Loaded geneset '{tag}': {len(genes)} genes from {gf.name}")

    if not genesets:
        print("[INFO] No valid gene lists loaded. Done.")
        return

    # Union genes across all contrasts (efficient: load expression ONCE)
    union_genes = []
    seen = set()
    for tag, glist in genesets.items():
        for g in glist:
            if g not in seen:
                seen.add(g)
                union_genes.append(g)

    print(f"[INFO] Union gene count across gene files: {len(union_genes)}")
    print("[INFO] Loading expression subset for UNION genes (merged cohort subjects only) ...")
    expr_union = load_expression_subset_subject_by_gene(
        expr_path=expr_path,
        subject_ids=ids_expr_mri,
        genes=union_genes,
        chunksize=args.chunksize
    )
    # expr_union columns are EXPR_<GENE> for genes that actually exist in expression file
    available_expr_cols = set(expr_union.columns)

    # Save union expression (optional but useful)
    union_out = out_dir / "selected_gene_expression_expr_mri_UNION_subject_by_gene.csv"
    expr_union.reset_index().to_csv(union_out, index=False)
    print("  Saved UNION expression:", union_out)

    # For each geneset, export selected expression + merged table
    for tag, glist in genesets.items():
        desired_cols = [f"EXPR_{g}" for g in glist]
        found_cols = [c for c in desired_cols if c in available_expr_cols]
        missing = [g for g in glist if f"EXPR_{g}" not in available_expr_cols]

        # Write found/missing lists
        (out_dir / f"gene_list_{tag}_found.txt").write_text("\n".join([c.replace("EXPR_", "") for c in found_cols]) + "\n")
        (out_dir / f"gene_list_{tag}_missing.txt").write_text("\n".join(missing) + "\n")

        print(f"[INFO] Geneset '{tag}': found={len(found_cols)} missing={len(missing)}")

        if len(found_cols) == 0:
            print(f"[WARN] No genes found in expression matrix for '{tag}'. Skipping exports.")
            continue

        expr_sub = expr_union[found_cols].copy()
        expr_sub_out = out_dir / f"selected_gene_expression_expr_mri_{tag}_subject_by_gene.csv"
        expr_sub.reset_index().to_csv(expr_sub_out, index=False)
        print("  Saved:", expr_sub_out)

        merged_with_expr = df_expr_mri.join(expr_sub, how="left")
        merged_out = out_dir / f"prepared_expr_mri_merged_with_selected_genes_{tag}.csv"
        merged_with_expr.reset_index().to_csv(merged_out, index=False)
        print("  Saved:", merged_out)

    print("[INFO] Done.")


if __name__ == "__main__":
    main()

