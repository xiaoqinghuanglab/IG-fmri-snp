"""Pairwise SNP-connectivity interaction analysis used in the paper workflow.

This module implements the manuscript-aligned imaging-genetics stage:
  1) load a union genotype matrix and combined SNP metadata
  2) load connectivity edges, covariates, and candidate resilience annotations
  3) run pairwise SNP x group interaction models for:
       - AsymAD vs TypicalAD
       - TypicalAD vs Control
       - AsymAD vs Control
  4) apply BH-FDR to interaction p-values
  5) export per-comparison full tables, significant hits, lead-SNP summaries,
     and lightweight diagnostics
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ALPHA = 0.05
MIN_TOTAL_N = 12
MIN_PER_GROUP_N = 3
MIN_CARRIERS_TOTAL = 2
NEAR_ZERO_VAR = 1e-12

COMPARISONS: List[Tuple[str, str, str]] = [
    ("AsymAD", "TypicalAD", "AsymAD_vs_TypAD"),
    ("TypicalAD", "Control", "TypAD_vs_Control"),
    ("AsymAD", "Control", "AsymAD_vs_Control"),
]

SNP_METADATA_REQUIRED_COLUMNS = ["SNP", "Case", "Control"]
SNP_METADATA_OPTIONAL_COLUMNS = ["CHR", "BP", "P", "OR"]
CANDIDATE_EDGE_COLUMNS = [
    "edge_id",
    "candidate_pattern_label",
    "strict_pattern_label",
    "adj_mean_Control",
    "adj_mean_AsymAD",
    "adj_mean_TypicalAD",
    "p_asym_vs_control",
    "p_typ_vs_control",
    "p_asym_vs_typ",
    "fdr_subtype_omnibus",
]


def normalize_subject_id(value) -> str | pd.NA:
    """Normalize ADNI subject identifiers to the shared 000_S_0000 style."""

    if pd.isna(value):
        return pd.NA

    text = str(value).strip()
    if text.startswith("sub-"):
        text = text[4:]

    match = re.match(r"^(\d{3})S(\d{4})$", text)
    if match:
        return f"{match.group(1)}_S_{match.group(2)}"

    match = re.match(r"^(\d{3})_S_(\d{4})$", text)
    if match:
        return text

    return text


def normalize_subtype(value) -> str | pd.NA:
    """Normalize subtype labels to Control / AsymAD / TypicalAD."""

    if pd.isna(value):
        return pd.NA

    text = str(value).strip()
    mapping = {
        "Control": "Control",
        "CN": "Control",
        "Normal": "Control",
        "Asym AD": "AsymAD",
        "AsymAD": "AsymAD",
        "Asym_AD": "AsymAD",
        "Typical AD": "TypicalAD",
        "TypicalAD": "TypicalAD",
        "Typical_AD": "TypicalAD",
        "TypAD": "TypicalAD",
    }
    return mapping.get(text, text)


def comparison_slug(case_group: str, control_group: str) -> str:
    """Build stable comparison labels while keeping TypAD filenames compact."""

    slug_map = {"TypicalAD": "TypAD"}
    case_slug = slug_map.get(case_group, case_group)
    control_slug = slug_map.get(control_group, control_group)
    return f"{case_slug}_vs_{control_slug}"


def safe_numeric_series(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def detect_subject_col(df: pd.DataFrame) -> str:
    for column in ["SubjectId", "Subject_ID", "subject_id", "subject", "ID", "id"]:
        if column in df.columns:
            return column
    return df.columns[0]


def reduce_to_full_rank(matrix: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    """Greedily keep columns only when they increase matrix rank."""

    keep_cols: List[str] = []
    dropped_cols: List[str] = []

    for column in matrix.columns:
        trial_cols = keep_cols + [column]
        trial = matrix[trial_cols].to_numpy(dtype=float)
        if np.linalg.matrix_rank(trial) > len(keep_cols):
            keep_cols.append(column)
        else:
            dropped_cols.append(column)

    return matrix[keep_cols].copy(), dropped_cols


def reduce_to_full_rank_prioritized(
    matrix: pd.DataFrame,
    priority_cols: Iterable[str],
) -> Tuple[pd.DataFrame, List[str]]:
    """Preserve key model columns whenever possible before dropping redundancy."""

    ordered_cols: List[str] = []
    seen = set()

    for column in priority_cols:
        if column in matrix.columns and column not in seen:
            ordered_cols.append(column)
            seen.add(column)

    for column in matrix.columns:
        if column not in seen:
            ordered_cols.append(column)
            seen.add(column)

    keep_cols: List[str] = []
    dropped_cols: List[str] = []

    for column in ordered_cols:
        trial_cols = keep_cols + [column]
        trial = matrix[trial_cols].to_numpy(dtype=float)
        if np.linalg.matrix_rank(trial) > len(keep_cols):
            keep_cols.append(column)
        else:
            dropped_cols.append(column)

    return matrix[keep_cols].copy(), dropped_cols


def add_fdr_columns(results_df: pd.DataFrame) -> pd.DataFrame:
    corrected = results_df.copy()

    for p_col, fdr_col in [
        ("p_interaction", "fdr_interaction"),
        ("p_snp_reference", "fdr_snp_reference"),
        ("p_snp_case", "fdr_snp_case"),
        ("p_group", "fdr_group"),
    ]:
        valid = corrected[p_col].notna()
        output = np.full(corrected.shape[0], np.nan)
        if valid.sum() > 0:
            output[valid] = multipletests(corrected.loc[valid, p_col], method="fdr_bh")[1]
        corrected[fdr_col] = output

    return corrected


def match_snp_columns(geno_columns: List[str], snp_rsids: List[str]) -> Tuple[List[Tuple[str, str]], List[str], List[str]]:
    rsid_to_col: Dict[str, str] = {}
    duplicated_rsids = set()

    for column in geno_columns:
        if column == "SubjectId":
            continue
        rsid = str(column).split("_")[0]
        if rsid not in rsid_to_col:
            rsid_to_col[rsid] = column
        else:
            duplicated_rsids.add(rsid)

    matched: List[Tuple[str, str]] = []
    missing: List[str] = []

    for rsid in snp_rsids:
        if rsid in rsid_to_col:
            matched.append((rsid, rsid_to_col[rsid]))
        else:
            missing.append(rsid)

    return matched, missing, sorted(duplicated_rsids)


def load_candidate_resilience_annotations(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Candidate resilience file not found: {path}")

    df = pd.read_csv(path)
    if "edge_id" not in df.columns:
        raise ValueError("Candidate resilience file must contain an edge_id column.")

    keep_cols = [column for column in CANDIDATE_EDGE_COLUMNS if column in df.columns]
    df = df[keep_cols].drop_duplicates(subset=["edge_id"]).copy()
    df["is_candidate_resilience"] = False
    if "candidate_pattern_label" in df.columns:
        df["is_candidate_resilience"] = df["candidate_pattern_label"].astype(str).str.startswith(
            "CandidateResilience_"
        )
    return df


def load_connectivity_wide(path: Path) -> Tuple[pd.DataFrame, List[str], str]:
    if not path.exists():
        raise FileNotFoundError(f"Connectivity file not found: {path}")

    conn = pd.read_csv(path)
    subject_col = detect_subject_col(conn)

    conn = conn.copy()
    conn["SubjectId"] = conn[subject_col].map(normalize_subject_id)

    edge_cols = [column for column in conn.columns if column not in [subject_col, "SubjectId"]]
    conn[edge_cols] = conn[edge_cols].apply(pd.to_numeric, errors="coerce")
    return conn, edge_cols, subject_col


def load_covariates(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Covariate file not found: {path}")

    cov = pd.read_csv(path)
    if "SubjectId" not in cov.columns:
        raise ValueError("Covariate file must contain SubjectId column.")

    cov = cov.copy()
    cov["SubjectId"] = cov["SubjectId"].map(normalize_subject_id)

    if "Subject_match_check" in cov.columns:
        cov = cov[cov["Subject_match_check"].fillna(True)].copy()

    required = ["SubjectId", "Age", "Sex", "Scan_type", "Manufacturer", "Subtype"]
    missing = [column for column in required if column not in cov.columns]
    if missing:
        raise ValueError(f"Covariate file missing required columns: {missing}")

    cov["Subtype_std"] = cov["Subtype"].map(normalize_subtype)
    cov["Age"] = safe_numeric_series(cov["Age"])
    cov["Sex_num"] = cov["Sex"].map({"F": 0.0, "M": 1.0})
    cov["Scan_type"] = cov["Scan_type"].fillna("Unknown").astype(str)
    cov["Manufacturer"] = cov["Manufacturer"].fillna("Unknown").astype(str)

    cov = cov.drop_duplicates(subset=["SubjectId"]).copy()
    return cov


def load_genotype_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Genotype file not found: {path}")

    df = pd.read_csv(path)
    subject_col = detect_subject_col(df)

    df = df.copy()
    df["SubjectId"] = df[subject_col].map(normalize_subject_id)

    snp_cols = [column for column in df.columns if column not in [subject_col, "SubjectId"]]
    df[snp_cols] = df[snp_cols].apply(pd.to_numeric, errors="coerce")
    return df[["SubjectId"] + snp_cols].drop_duplicates(subset=["SubjectId"]).copy()


def load_snp_metadata(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"SNP metadata file not found: {path}")

    df = pd.read_csv(path)
    missing_required = [column for column in SNP_METADATA_REQUIRED_COLUMNS if column not in df.columns]
    if missing_required:
        raise ValueError(f"{path.name} missing required columns: {missing_required}")

    df = df.copy()
    df["SNP"] = df["SNP"].astype(str).str.strip()
    df["Case_std"] = df["Case"].map(normalize_subtype)
    df["Control_std"] = df["Control"].map(normalize_subtype)
    df["source_comparison"] = df.apply(
        lambda row: comparison_slug(str(row["Case_std"]), str(row["Control_std"])),
        axis=1,
    )

    for column in SNP_METADATA_OPTIONAL_COLUMNS:
        if column in df.columns:
            df[column] = safe_numeric_series(df[column])

    source_summary = (
        df.groupby("SNP")["source_comparison"]
        .agg(lambda values: ";".join(sorted(set(map(str, values)))))
        .reset_index()
        .rename(columns={"source_comparison": "source_comparisons"})
    )

    n_source_summary = (
        df.groupby("SNP")["source_comparison"]
        .nunique()
        .reset_index()
        .rename(columns={"source_comparison": "n_source_comparisons"})
    )

    sort_cols = ["SNP"]
    if "P" in df.columns:
        sort_cols = ["P", "SNP"]

    meta_best = df.sort_values(sort_cols, na_position="last").drop_duplicates(subset=["SNP"], keep="first").copy()
    meta_best = meta_best.merge(source_summary, on="SNP", how="left")
    meta_best = meta_best.merge(n_source_summary, on="SNP", how="left")
    return meta_best


def prepare_master_analysis_df(
    connectivity: pd.DataFrame,
    covariates: pd.DataFrame,
    genotypes: pd.DataFrame,
) -> pd.DataFrame:
    master = connectivity.merge(covariates, on="SubjectId", how="inner")
    master = master.merge(genotypes, on="SubjectId", how="left")
    return master


def prepare_comparison_df(
    master: pd.DataFrame,
    case_group: str,
    reference_group: str,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """Filter to one pairwise contrast and build nuisance covariates."""

    df = master.copy()
    df = df[df["Subtype_std"].isin([case_group, reference_group])].copy()
    df["Group"] = (df["Subtype_std"] == case_group).astype(float)
    df = df.dropna(subset=["Age", "Sex_num"]).copy()

    scan_dummies = pd.get_dummies(df["Scan_type"], prefix="Scan", drop_first=True, dtype=float)
    manufacturer_dummies = pd.get_dummies(df["Manufacturer"], prefix="Manufacturer", drop_first=True, dtype=float)
    dummy_df = pd.concat([scan_dummies, manufacturer_dummies], axis=1)

    if dummy_df.shape[1] > 0:
        dummy_df, dropped_dummy_cols = reduce_to_full_rank(dummy_df)
        df = pd.concat([df.reset_index(drop=True), dummy_df.reset_index(drop=True)], axis=1)
    else:
        dropped_dummy_cols = []
        df = df.reset_index(drop=True)

    nuisance_cols = ["Age", "Sex_num"] + dummy_df.columns.tolist()
    return df, nuisance_cols, dropped_dummy_cols


def fit_one_snp_all_edges(
    df_comp: pd.DataFrame,
    edge_cols: List[str],
    nuisance_cols: List[str],
    snp_col: str,
    case_group: str,
    reference_group: str,
) -> Tuple[Optional[dict], Optional[dict]]:
    """Fit one SNP across all edges for a single two-group comparison."""

    df = df_comp.copy()
    df["SNP_raw"] = safe_numeric_series(df[snp_col])
    df = df.loc[df["SNP_raw"].notna()].copy()

    n_total = df.shape[0]
    n_case = int((df["Subtype_std"] == case_group).sum())
    n_reference = int((df["Subtype_std"] == reference_group).sum())

    if n_total < MIN_TOTAL_N:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": f"TooFewSubjects_total<{MIN_TOTAL_N}",
        }

    if min(n_case, n_reference) < MIN_PER_GROUP_N:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": f"TooFewSubjectsInOneGroup_min<{MIN_PER_GROUP_N}",
        }

    case_vals = df.loc[df["Subtype_std"] == case_group, "SNP_raw"].dropna()
    reference_vals = df.loc[df["Subtype_std"] == reference_group, "SNP_raw"].dropna()

    if case_vals.nunique(dropna=True) < 2:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": "NoWithinCaseVariation",
        }

    if reference_vals.nunique(dropna=True) < 2:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": "NoWithinReferenceVariation",
        }

    n_carriers_total = int((df["SNP_raw"] > 0).sum())
    n_carriers_case = int((case_vals > 0).sum())
    n_carriers_reference = int((reference_vals > 0).sum())

    if n_carriers_total < MIN_CARRIERS_TOTAL:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": f"TooFewCarriersTotal<{MIN_CARRIERS_TOTAL}",
        }

    if float(np.nanvar(df["SNP_raw"].to_numpy(dtype=float), ddof=1)) <= NEAR_ZERO_VAR:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": "NearZeroSNPVariance",
        }

    df["SNP_c"] = df["SNP_raw"] - df["SNP_raw"].mean()
    df["SNP_x_Group"] = df["SNP_c"] * df["Group"]

    design = pd.concat(
        [
            pd.Series(1.0, index=df.index, name="Intercept"),
            df[["SNP_c", "Group", "SNP_x_Group"]],
            df[nuisance_cols],
        ],
        axis=1,
    )

    design, dropped_full_cols = reduce_to_full_rank_prioritized(
        design,
        priority_cols=["Intercept", "SNP_c", "Group", "SNP_x_Group", "Age", "Sex_num"],
    )

    for required_col in ["Intercept", "SNP_c", "Group", "SNP_x_Group"]:
        if required_col not in design.columns:
            return None, {
                "snp_col": snp_col,
                "n_total": n_total,
                "n_case": n_case,
                "n_reference": n_reference,
                "reason": f"DroppedRequiredColumn_{required_col}",
            }

    design_values = design.to_numpy(dtype=float)
    rank = np.linalg.matrix_rank(design_values)
    p = design.shape[1]
    df_resid = n_total - rank

    if rank < p:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": "RankDeficientDesign",
        }

    if df_resid <= 0:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": "NonPositiveResidualDF",
        }

    response = df[edge_cols].apply(pd.to_numeric, errors="coerce")
    usable_edges = response.notna().all(axis=0)
    usable_edge_names = usable_edges[usable_edges].index.tolist()
    if not usable_edge_names:
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": "AllEdgesMissingForSNPSample",
        }

    response = response[usable_edge_names].copy()
    response_values = response.to_numpy(dtype=float)

    try:
        xtx = design_values.T @ design_values
        xtx_inv = np.linalg.pinv(xtx)
        betas = xtx_inv @ design_values.T @ response_values
        y_hat = design_values @ betas
        residuals = response_values - y_hat

        sse = np.sum(residuals ** 2, axis=0)
        sigma2 = sse / df_resid
        y_mean = np.mean(response_values, axis=0)
        sst = np.sum((response_values - y_mean) ** 2, axis=0)

        with np.errstate(divide="ignore", invalid="ignore"):
            r2 = 1.0 - (sse / sst)
        r2 = np.where(np.isfinite(r2), r2, np.nan)
        adj_r2 = 1.0 - (1.0 - r2) * ((n_total - 1.0) / df_resid)

        col_idx = {column: i for i, column in enumerate(design.columns)}

        def coef_stats(colname: str):
            idx = col_idx[colname]
            beta = betas[idx, :]
            se = np.sqrt(np.maximum(xtx_inv[idx, idx] * sigma2, 0))
            tval = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
            pval = 2 * stats.t.sf(np.abs(tval), df_resid)
            return beta, se, tval, pval

        def combo_stats(weights: np.ndarray):
            vector = np.asarray(weights, dtype=float).reshape(-1, 1)
            beta = (vector.T @ betas).ravel()
            variance = float((vector.T @ xtx_inv @ vector).item())
            se = np.sqrt(np.maximum(variance * sigma2, 0))
            tval = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
            pval = 2 * stats.t.sf(np.abs(tval), df_resid)
            return beta, se, tval, pval

        beta_snp_ref, se_snp_ref, t_snp_ref, p_snp_ref = coef_stats("SNP_c")
        beta_group, se_group, t_group, p_group = coef_stats("Group")
        beta_interaction, se_interaction, t_interaction, p_interaction = coef_stats("SNP_x_Group")

        combo = np.zeros(design.shape[1], dtype=float)
        combo[col_idx["SNP_c"]] = 1.0
        combo[col_idx["SNP_x_Group"]] = 1.0
        beta_snp_case, se_snp_case, t_snp_case, p_snp_case = combo_stats(combo)

        mean_edge_case = df.loc[df["Subtype_std"] == case_group, usable_edge_names].mean(axis=0)
        mean_edge_reference = df.loc[df["Subtype_std"] == reference_group, usable_edge_names].mean(axis=0)

        return {
            "usable_edges": usable_edge_names,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "n_carriers_total": n_carriers_total,
            "n_carriers_case": n_carriers_case,
            "n_carriers_reference": n_carriers_reference,
            "df_resid": df_resid,
            "rank": rank,
            "dropped_full_cols": ";".join(dropped_full_cols) if dropped_full_cols else "",
            "beta_snp_reference": beta_snp_ref,
            "se_snp_reference": se_snp_ref,
            "t_snp_reference": t_snp_ref,
            "p_snp_reference": p_snp_ref,
            "beta_group": beta_group,
            "se_group": se_group,
            "t_group": t_group,
            "p_group": p_group,
            "beta_interaction": beta_interaction,
            "se_interaction": se_interaction,
            "t_interaction": t_interaction,
            "p_interaction": p_interaction,
            "beta_snp_case": beta_snp_case,
            "se_snp_case": se_snp_case,
            "t_snp_case": t_snp_case,
            "p_snp_case": p_snp_case,
            "r2": r2,
            "adj_r2": adj_r2,
            "mean_edge_case": mean_edge_case,
            "mean_edge_reference": mean_edge_reference,
        }, None
    except Exception as exc:  # pragma: no cover - protective batch-run guard
        return None, {
            "snp_col": snp_col,
            "n_total": n_total,
            "n_case": n_case,
            "n_reference": n_reference,
            "reason": f"ModelError:{exc}",
        }


def run_pairwise_analysis(
    connectivity_path: Path,
    covariates_path: Path,
    genetic_path: Path,
    snp_metadata_path: Path,
    candidate_edges_path: Path,
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n======================")
    print("PAIRWISE SNP-CONNECTIVITY ANALYSIS")
    print("======================")
    print(f"Connectivity file        : {connectivity_path}")
    print(f"Covariate file           : {covariates_path}")
    print(f"Genotype matrix          : {genetic_path}")
    print(f"SNP metadata             : {snp_metadata_path}")
    print(f"Candidate resilience     : {candidate_edges_path}")
    print(f"Output dir               : {output_dir}")

    resilience_annot = load_candidate_resilience_annotations(candidate_edges_path)
    connectivity, all_edge_cols, _ = load_connectivity_wide(connectivity_path)
    covariates = load_covariates(covariates_path)
    genotype_matrix = load_genotype_matrix(genetic_path)
    snp_metadata = load_snp_metadata(snp_metadata_path)

    print("\n[Input shapes]")
    print(f"Connectivity raw shape   : {connectivity.shape}")
    print(f"Covariate raw shape      : {covariates.shape}")
    print(f"Genotype raw shape       : {genotype_matrix.shape}")
    print(f"SNP metadata rows        : {snp_metadata.shape[0]}")
    print(f"Connectivity edge count  : {len(all_edge_cols)}")
    print(f"Candidate resilience edge count: {resilience_annot['edge_id'].nunique()}")

    matched_snps, missing_snps, duplicated_rsids = match_snp_columns(
        genotype_matrix.columns.tolist(),
        snp_metadata["SNP"].tolist(),
    )
    matched_snp_df = pd.DataFrame(matched_snps, columns=["SNP", "Genotype_Column"])
    snp_metadata = snp_metadata.merge(matched_snp_df, on="SNP", how="inner").copy()

    print("\n[Union genotype / SNP metadata]")
    print(f"Matched SNPs             : {snp_metadata.shape[0]}")
    print(f"Missing SNPs             : {len(missing_snps)}")
    print(f"Duplicated rsIDs         : {len(duplicated_rsids)}")

    master = prepare_master_analysis_df(connectivity, covariates, genotype_matrix)

    print("\n[Merged master data]")
    print(f"Master shape             : {master.shape}")
    print(master["Subtype_std"].value_counts(dropna=False).to_string())

    for case_group, reference_group, comparison_name in COMPARISONS:
        comp_dir = output_dir / comparison_name
        comp_dir.mkdir(parents=True, exist_ok=True)

        print("\n" + "=" * 80)
        print(f"RUNNING COMPARISON: {comparison_name}")
        print("=" * 80)

        df_comp, nuisance_cols, dropped_nuisance_cols = prepare_comparison_df(master, case_group, reference_group)

        print("\n[Comparison base data]")
        print(f"Shape                   : {df_comp.shape}")
        print(f"All connectivity edges  : {len(all_edge_cols)}")
        print(df_comp["Subtype_std"].value_counts(dropna=False).to_string())

        if dropped_nuisance_cols:
            print("[Dropped redundant nuisance columns]")
            for column in dropped_nuisance_cols:
                print(f"  - {column}")

        if df_comp.shape[0] < MIN_TOTAL_N:
            print("[SKIP] Too few subjects overall for this comparison.")
            continue

        if df_comp["Subtype_std"].nunique() < 2:
            print("[SKIP] Comparison does not contain both groups.")
            continue

        if snp_metadata.empty:
            print("[SKIP] No SNPs matched genotype matrix columns.")
            continue

        results_chunks: List[pd.DataFrame] = []
        skipped_snps: List[dict] = []
        total_snps = snp_metadata.shape[0]

        for i, row in snp_metadata.reset_index(drop=True).iterrows():
            snp = row["SNP"]
            snp_col = row["Genotype_Column"]

            if snp_col not in df_comp.columns:
                skipped_snps.append(
                    {
                        "SNP": snp,
                        "genotype_column": snp_col,
                        "reason": "GenotypeColumnMissingInMasterTable",
                    }
                )
                continue

            fit, skip = fit_one_snp_all_edges(
                df_comp=df_comp,
                edge_cols=all_edge_cols,
                nuisance_cols=nuisance_cols,
                snp_col=snp_col,
                case_group=case_group,
                reference_group=reference_group,
            )

            if fit is None:
                out = {"SNP": snp, "genotype_column": snp_col}
                if skip is not None:
                    out.update(skip)
                skipped_snps.append(out)
            else:
                usable_edges = fit["usable_edges"]
                result = pd.DataFrame({"edge_id": usable_edges})
                result["comparison"] = comparison_name
                result["case_group"] = case_group
                result["reference_group"] = reference_group
                result["SNP"] = snp
                result["genotype_column"] = snp_col

                for column in ["CHR", "BP", "P", "OR", "Case_std", "Control_std", "source_comparison"]:
                    if column in row.index:
                        out_col = {
                            "P": "GWAS_P",
                            "OR": "GWAS_OR",
                            "Case_std": "GWAS_Case",
                            "Control_std": "GWAS_Control",
                        }.get(column, column)
                        result[out_col] = row[column]

                result["source_comparisons"] = row.get("source_comparisons", np.nan)
                result["n_source_comparisons"] = row.get("n_source_comparisons", np.nan)
                result["n_subjects"] = fit["n_total"]
                result["n_case"] = fit["n_case"]
                result["n_reference"] = fit["n_reference"]
                result["n_carriers_total"] = fit["n_carriers_total"]
                result["n_carriers_case"] = fit["n_carriers_case"]
                result["n_carriers_reference"] = fit["n_carriers_reference"]
                result["df_resid"] = fit["df_resid"]
                result["rank"] = fit["rank"]
                result["dropped_full_cols"] = fit["dropped_full_cols"]
                result["beta_snp_reference"] = fit["beta_snp_reference"]
                result["se_snp_reference"] = fit["se_snp_reference"]
                result["t_snp_reference"] = fit["t_snp_reference"]
                result["p_snp_reference"] = fit["p_snp_reference"]
                result["beta_group"] = fit["beta_group"]
                result["se_group"] = fit["se_group"]
                result["t_group"] = fit["t_group"]
                result["p_group"] = fit["p_group"]
                result["beta_interaction"] = fit["beta_interaction"]
                result["se_interaction"] = fit["se_interaction"]
                result["t_interaction"] = fit["t_interaction"]
                result["p_interaction"] = fit["p_interaction"]
                result["beta_snp_case"] = fit["beta_snp_case"]
                result["se_snp_case"] = fit["se_snp_case"]
                result["t_snp_case"] = fit["t_snp_case"]
                result["p_snp_case"] = fit["p_snp_case"]
                result["r2"] = fit["r2"]
                result["adj_r2"] = fit["adj_r2"]
                result["mean_edge_case"] = result["edge_id"].map(fit["mean_edge_case"].to_dict())
                result["mean_edge_reference"] = result["edge_id"].map(fit["mean_edge_reference"].to_dict())
                results_chunks.append(result)

            if ((i + 1) % 200 == 0) or ((i + 1) == total_snps):
                print(f"  Processed SNPs: {i + 1}/{total_snps}")

        skipped_df = pd.DataFrame(skipped_snps)
        if not results_chunks:
            print("[SKIP] No valid SNP-edge models were fit for this comparison.")
            continue

        results_df = pd.concat(results_chunks, ignore_index=True)
        results_df = add_fdr_columns(results_df)
        results_df = results_df.merge(resilience_annot, on="edge_id", how="left")
        results_df["is_candidate_resilience"] = results_df["is_candidate_resilience"].fillna(False)
        results_df = results_df.sort_values(
            ["fdr_interaction", "p_interaction", "p_snp_case", "p_snp_reference"],
            na_position="last",
        ).reset_index(drop=True)

        full_out = comp_dir / f"{comparison_name}_full_results.csv.gz"

        results_df.to_csv(full_out, index=False, compression="gzip")

        print(f"\n[Done] {comparison_name}")
        print(f"  Full results : {full_out}")
        print(f"  Result rows  : {results_df.shape[0]}")

    print("\nAll comparisons finished.")
