#!/usr/bin/env python3

from pathlib import Path
import re
import warnings
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# =========================================================
# DEFAULT PATHS
# =========================================================
BASE_DIR = Path("~/IG-fmri-snp")
DEFAULT_CONNECTIVITY_FILE = BASE_DIR / "msdl_all_subjects_connectivity_edges.csv"
DEFAULT_COVARIATE_FILE = BASE_DIR / "covariate_file.csv"
DEFAULT_OUTPUT_DIR = BASE_DIR / "Model2_MSDL_Phenotype"

# =========================================================
# CONFIG
# =========================================================
ALPHA = 0.05
MIN_TOTAL_N = 20
MIN_PER_GROUP_N = 3
NEAR_ZERO_VAR = 1e-12

# Candidate-screen thresholds
CANDIDATE_ALPHA = 0.05
CANDIDATE_MIN_ABS_BETA = 0.0  # keep 0.0 unless you want an effect-size floor

# =========================================================
# HELPERS
# =========================================================
def normalize_subject_id(x):
    """
    Convert:
      002_S_1155 -> 002_S_1155
      sub-002S1155 -> 002_S_1155
      002S1155 -> 002_S_1155
    """
    if pd.isna(x):
        return pd.NA

    s = str(x).strip()

    if s.startswith("sub-"):
        s = s[4:]

    m = re.match(r"^(\d{3})S(\d{4})$", s)
    if m:
        return f"{m.group(1)}_S_{m.group(2)}"

    m = re.match(r"^(\d{3})_S_(\d{4})$", s)
    if m:
        return s

    return s


def normalize_subtype(x):
    """
    Standardize subtype labels to:
      Control, AsymAD, TypicalAD
    """
    if pd.isna(x):
        return pd.NA

    s = str(x).strip()

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

    return mapping.get(s, s)


def safe_numeric_series(s):
    return pd.to_numeric(s, errors="coerce")


def print_counts_block(title, s):
    print(f"\n{title}")
    print("-" * len(title))
    if len(s) == 0:
        print("  (none)")
    else:
        print(s.to_string())


def reduce_to_full_rank(X):
    """
    Greedily keep columns only if they increase matrix rank.
    Returns:
        X_reduced, dropped_columns
    """
    keep_cols = []
    dropped_cols = []

    for col in X.columns:
        trial_cols = keep_cols + [col]
        trial = X[trial_cols].to_numpy(dtype=float)
        if np.linalg.matrix_rank(trial) > len(keep_cols):
            keep_cols.append(col)
        else:
            dropped_cols.append(col)

    return X[keep_cols].copy(), dropped_cols


def joint_f_test(beta, XtX_inv, sigma2, idx_a, idx_b, df_resid):
    """
    Joint test H0: beta[idx_a] = 0 and beta[idx_b] = 0
    """
    p = len(beta)
    R = np.zeros((2, p), dtype=float)
    R[0, idx_a] = 1.0
    R[1, idx_b] = 1.0

    Rb = R @ beta
    middle = R @ XtX_inv @ R.T

    try:
        middle_inv = np.linalg.pinv(middle)
        q = 2
        F = float((Rb.T @ middle_inv @ Rb) / (q * sigma2))
        if not np.isfinite(F):
            return np.nan
        return float(stats.f.sf(F, q, df_resid))
    except Exception:
        return np.nan


def linear_combo_stats(beta, XtX_inv, sigma2, weights, df_resid):
    """
    For a contrast c'beta, return estimate, se, t, p
    """
    w = np.asarray(weights, dtype=float).reshape(-1, 1)
    est = float((w.T @ beta.reshape(-1, 1)).item())
    var = float((w.T @ XtX_inv @ w).item()) * sigma2
    var = max(var, 0.0)
    se = np.sqrt(var)
    if se <= 0:
        return est, np.nan, np.nan, np.nan
    tval = est / se
    pval = 2 * stats.t.sf(np.abs(tval), df_resid)
    return est, se, tval, pval


def build_base_design(merged):
    """
    Build model:
      edge ~ C(Subtype) + Age + Sex_num + C(Scan_type) + C(Manufacturer)

    Reference subtype = Control
    """
    df = merged.copy()

    df["Age"] = safe_numeric_series(df["Age"])
    df["Sex_num"] = df["Sex"].map({"F": 0.0, "M": 1.0})
    df["Scan_type"] = df["Scan_type"].fillna("Unknown").astype(str)
    df["Manufacturer"] = df["Manufacturer"].fillna("Unknown").astype(str)
    df["Subtype_std"] = df["Subtype"].map(normalize_subtype)

    df = df[df["Subtype_std"].isin(["Control", "AsymAD", "TypicalAD"])].copy()
    df = df.dropna(subset=["Age", "Sex_num", "Subtype_std"]).copy()

    df["Subtype_AsymAD"] = (df["Subtype_std"] == "AsymAD").astype(float)
    df["Subtype_TypicalAD"] = (df["Subtype_std"] == "TypicalAD").astype(float)

    scan_dummies = pd.get_dummies(df["Scan_type"], prefix="Scan", drop_first=True, dtype=float)
    manuf_dummies = pd.get_dummies(df["Manufacturer"], prefix="Manufacturer", drop_first=True, dtype=float)

    X = pd.concat(
        [
            pd.Series(1.0, index=df.index, name="Intercept"),
            df[["Subtype_AsymAD", "Subtype_TypicalAD", "Age", "Sex_num"]],
            scan_dummies,
            manuf_dummies,
        ],
        axis=1,
    )

    return df, X


def fit_edge_ols(X, y, subtype_series):
    """
    Fit OLS for one edge.
    Returns dict or a reason to skip.
    """
    mask = y.notna()
    X_sub = X.loc[mask].copy()
    y_sub = safe_numeric_series(y.loc[mask]).astype(float)
    subtype_sub = subtype_series.loc[mask].copy()

    n_total = int(mask.sum())
    n_control = int((subtype_sub == "Control").sum())
    n_asym = int((subtype_sub == "AsymAD").sum())
    n_typ = int((subtype_sub == "TypicalAD").sum())

    if n_total < MIN_TOTAL_N:
        return None, {
            "n_total": n_total,
            "n_Control": n_control,
            "n_AsymAD": n_asym,
            "n_TypicalAD": n_typ,
            "reason": f"TooFewSubjects_total<{MIN_TOTAL_N}",
        }

    if min(n_control, n_asym, n_typ) < MIN_PER_GROUP_N:
        return None, {
            "n_total": n_total,
            "n_Control": n_control,
            "n_AsymAD": n_asym,
            "n_TypicalAD": n_typ,
            "reason": f"TooFewSubjectsInOneGroup_min<{MIN_PER_GROUP_N}",
        }

    y_std = float(np.nanstd(y_sub.to_numpy(dtype=float), ddof=1))
    if (not np.isfinite(y_std)) or y_std <= NEAR_ZERO_VAR:
        return None, {
            "n_total": n_total,
            "n_Control": n_control,
            "n_AsymAD": n_asym,
            "n_TypicalAD": n_typ,
            "reason": "NearZeroVarianceEdge",
        }

    X_sub, dropped_edge_cols = reduce_to_full_rank(X_sub)
    Xv = X_sub.to_numpy(dtype=float)
    yv = y_sub.to_numpy(dtype=float)

    rank = np.linalg.matrix_rank(Xv)
    p = Xv.shape[1]
    df_resid = n_total - rank

    if rank < p or df_resid <= 0:
        return None, {
            "n_total": n_total,
            "n_Control": n_control,
            "n_AsymAD": n_asym,
            "n_TypicalAD": n_typ,
            "reason": "RankDeficientDesign",
        }

    try:
        XtX = Xv.T @ Xv
        XtX_inv = np.linalg.pinv(XtX)
        beta = XtX_inv @ Xv.T @ yv
        yhat = Xv @ beta
        resid = yv - yhat

        sse = float(np.sum(resid ** 2))
        sigma2 = sse / df_resid

        y_mean = float(np.mean(yv))
        sst = float(np.sum((yv - y_mean) ** 2))
        r2 = np.nan
        adj_r2 = np.nan
        if sst > 0:
            r2 = 1.0 - (sse / sst)
            adj_r2 = 1.0 - (1.0 - r2) * ((n_total - 1.0) / df_resid)

        col_index = {c: i for i, c in enumerate(X_sub.columns)}

        if ("Subtype_AsymAD" not in col_index) or ("Subtype_TypicalAD" not in col_index):
            return None, {
                "n_total": n_total,
                "n_Control": n_control,
                "n_AsymAD": n_asym,
                "n_TypicalAD": n_typ,
                "reason": "DroppedSubtypeColumnDueToCollinearity",
            }

        idx_asym = col_index["Subtype_AsymAD"]
        idx_typ = col_index["Subtype_TypicalAD"]

        se_beta = np.sqrt(np.maximum(np.diag(XtX_inv) * sigma2, 0.0))

        beta_asym = float(beta[idx_asym])
        se_asym = float(se_beta[idx_asym])
        t_asym = beta_asym / se_asym if se_asym > 0 else np.nan
        p_asym = 2 * stats.t.sf(np.abs(t_asym), df_resid) if np.isfinite(t_asym) else np.nan

        beta_typ = float(beta[idx_typ])
        se_typ = float(se_beta[idx_typ])
        t_typ = beta_typ / se_typ if se_typ > 0 else np.nan
        p_typ = 2 * stats.t.sf(np.abs(t_typ), df_resid) if np.isfinite(t_typ) else np.nan

        w = np.zeros(len(beta), dtype=float)
        w[idx_asym] = 1.0
        w[idx_typ] = -1.0
        beta_asym_vs_typ, se_asym_vs_typ, t_asym_vs_typ, p_asym_vs_typ = linear_combo_stats(
            beta, XtX_inv, sigma2, w, df_resid
        )

        p_subtype_omnibus = joint_f_test(beta, XtX_inv, sigma2, idx_asym, idx_typ, df_resid)

        x_mean = X_sub.mean(axis=0).to_dict()
        x_control = {k: float(v) for k, v in x_mean.items()}
        x_asym_mean = {k: float(v) for k, v in x_mean.items()}
        x_typ_mean = {k: float(v) for k, v in x_mean.items()}

        x_control["Intercept"] = 1.0
        x_asym_mean["Intercept"] = 1.0
        x_typ_mean["Intercept"] = 1.0

        x_control["Subtype_AsymAD"] = 0.0
        x_control["Subtype_TypicalAD"] = 0.0

        x_asym_mean["Subtype_AsymAD"] = 1.0
        x_asym_mean["Subtype_TypicalAD"] = 0.0

        x_typ_mean["Subtype_AsymAD"] = 0.0
        x_typ_mean["Subtype_TypicalAD"] = 1.0

        col_order = X_sub.columns.tolist()
        b_ser = pd.Series(beta, index=col_order)

        adj_mean_control = float(sum(x_control[c] * b_ser[c] for c in col_order))
        adj_mean_asym = float(sum(x_asym_mean[c] * b_ser[c] for c in col_order))
        adj_mean_typ = float(sum(x_typ_mean[c] * b_ser[c] for c in col_order))

        out = {
            "n_total": n_total,
            "n_Control": n_control,
            "n_AsymAD": n_asym,
            "n_TypicalAD": n_typ,
            "n_covariates_used": X_sub.shape[1],
            "df_resid": df_resid,
            "rank": rank,
            "r2": r2,
            "adj_r2": adj_r2,
            "beta_asym_vs_control": beta_asym,
            "se_asym_vs_control": se_asym,
            "t_asym_vs_control": t_asym,
            "p_asym_vs_control": p_asym,
            "beta_typ_vs_control": beta_typ,
            "se_typ_vs_control": se_typ,
            "t_typ_vs_control": t_typ,
            "p_typ_vs_control": p_typ,
            "beta_asym_vs_typ": beta_asym_vs_typ,
            "se_asym_vs_typ": se_asym_vs_typ,
            "t_asym_vs_typ": t_asym_vs_typ,
            "p_asym_vs_typ": p_asym_vs_typ,
            "p_subtype_omnibus": p_subtype_omnibus,
            "adj_mean_Control": adj_mean_control,
            "adj_mean_AsymAD": adj_mean_asym,
            "adj_mean_TypicalAD": adj_mean_typ,
            "edge_sd": y_std,
            "dropped_edge_design_cols": ";".join(dropped_edge_cols) if dropped_edge_cols else "",
        }
        return out, None

    except Exception as e:
        return None, {
            "n_total": n_total,
            "n_Control": n_control,
            "n_AsymAD": n_asym,
            "n_TypicalAD": n_typ,
            "reason": f"ModelError:{str(e)}",
        }


def add_fdr_columns(df):
    df = df.copy()

    for pcol, fdrcol in [
        ("p_subtype_omnibus", "fdr_subtype_omnibus"),
        ("p_asym_vs_control", "fdr_asym_vs_control"),
        ("p_typ_vs_control", "fdr_typ_vs_control"),
        ("p_asym_vs_typ", "fdr_asym_vs_typ"),
    ]:
        valid = df[pcol].notna()
        out = np.full(df.shape[0], np.nan)
        if valid.sum() > 0:
            out[valid] = multipletests(df.loc[valid, pcol], method="fdr_bh")[1]
        df[fdrcol] = out

    return df


def classify_pattern_strict(row, alpha=ALPHA):
    """
    Strict pattern logic based on FDR thresholds + omnibus signal.
    """
    m_c = row["adj_mean_Control"]
    m_a = row["adj_mean_AsymAD"]
    m_t = row["adj_mean_TypicalAD"]

    fa = row["fdr_asym_vs_control"]
    ft = row["fdr_typ_vs_control"]
    fat = row["fdr_asym_vs_typ"]
    fo = row["fdr_subtype_omnibus"]

    if pd.isna(fo) or fo >= alpha:
        return "MixedOrUnclear_NoOmnibusSignal"

    sig_a = pd.notna(fa) and (fa < alpha)
    sig_t = pd.notna(ft) and (ft < alpha)
    sig_at = pd.notna(fat) and (fat < alpha)

    if (not sig_a) and sig_t and sig_at:
        if (m_a >= m_t) and (m_c >= m_t):
            return "Resilience_AsymAD≈Control>TypicalAD"
        if (m_a <= m_t) and (m_c <= m_t):
            return "Resilience_AsymAD≈Control<TypicalAD"
        return "Resilience_AsymAD≈Control_DirectionMixed"

    if sig_a and sig_t and (not sig_at):
        if (m_c >= m_a) and (m_c >= m_t):
            return "NonResilience_Control>AsymAD≈TypicalAD"
        if (m_c <= m_a) and (m_c <= m_t):
            return "NonResilience_Control<AsymAD≈TypicalAD"
        return "NonResilience_ControlDifferent_DirectionMixed"

    if sig_a and sig_at:
        if (m_a > m_c) and (m_a > m_t):
            return "Compensation_AsymADHighest"
        if (m_a < m_c) and (m_a < m_t):
            return "Compensation_AsymADLowest"
        return "Compensation_AsymADDistinct"

    means = {"Control": m_c, "AsymAD": m_a, "TypicalAD": m_t}
    ordered = sorted(means.items(), key=lambda x: x[1], reverse=True)
    order_str = ">".join([k for k, _ in ordered])
    return f"MixedOrUnclear_{order_str}"


def classify_pattern_candidate(row, alpha=CANDIDATE_ALPHA, min_abs_beta=CANDIDATE_MIN_ABS_BETA):
    """
    Candidate pattern logic using RAW p-values + adjusted mean ordering.
    This is for screening edges for downstream resilience-focused SNP analysis.
    """
    m_c = row["adj_mean_Control"]
    m_a = row["adj_mean_AsymAD"]
    m_t = row["adj_mean_TypicalAD"]

    pa = row["p_asym_vs_control"]
    pt = row["p_typ_vs_control"]
    pat = row["p_asym_vs_typ"]

    ba = row["beta_asym_vs_control"]
    bt = row["beta_typ_vs_control"]
    bat = row["beta_asym_vs_typ"]

    sig_a = pd.notna(pa) and (pa < alpha) and (abs(ba) >= min_abs_beta)
    sig_t = pd.notna(pt) and (pt < alpha) and (abs(bt) >= min_abs_beta)
    sig_at = pd.notna(pat) and (pat < alpha) and (abs(bat) >= min_abs_beta)

    # Candidate resilience: Asym ~ Control, both different from Typical
    if (not sig_a) and sig_t and sig_at:
        if (m_a >= m_t) and (m_c >= m_t):
            return "CandidateResilience_AsymAD≈Control>TypicalAD"
        if (m_a <= m_t) and (m_c <= m_t):
            return "CandidateResilience_AsymAD≈Control<TypicalAD"
        return "CandidateResilience_AsymAD≈Control_DirectionMixed"

    # Candidate non-resilience/pathology-like: Asym ~ Typical, both different from Control
    if sig_a and sig_t and (not sig_at):
        if (m_c >= m_a) and (m_c >= m_t):
            return "CandidateNonResilience_Control>AsymAD≈TypicalAD"
        if (m_c <= m_a) and (m_c <= m_t):
            return "CandidateNonResilience_Control<AsymAD≈TypicalAD"
        return "CandidateNonResilience_ControlDifferent_DirectionMixed"

    # Candidate compensation: Asym differs from both
    if sig_a and sig_at:
        if (m_a > m_c) and (m_a > m_t):
            return "CandidateCompensation_AsymADHighest"
        if (m_a < m_c) and (m_a < m_t):
            return "CandidateCompensation_AsymADLowest"
        return "CandidateCompensation_AsymADDistinct"

    means = {"Control": m_c, "AsymAD": m_a, "TypicalAD": m_t}
    ordered = sorted(means.items(), key=lambda x: x[1], reverse=True)
    order_str = ">".join([k for k, _ in ordered])
    return f"CandidateMixedOrUnclear_{order_str}"


def write_diagnostic_report(
    out_path,
    connectivity_file,
    covariate_file,
    output_dir,
    conn_shape,
    cov_shape,
    merged_before_shape,
    merged_after_shape,
    analysis_shape,
    n_edges,
    subtype_counts_before,
    subtype_counts_after,
    X_base_shape,
    X_base_rank,
    dropped_rank_cols,
    edge_n_summary,
    near_zero_var_edges,
    warnings_list,
    results_df,
    skipped_df,
):
    lines = []
    lines.append("MODEL 2: MSDL phenotype-only subtype model\n")
    lines.append(f"Connectivity file: {connectivity_file}")
    lines.append(f"Covariate file: {covariate_file}")
    lines.append(f"Output dir: {output_dir}\n")

    lines.append("[Input shapes]")
    lines.append(f"Connectivity raw shape: {conn_shape}")
    lines.append(f"Covariate raw shape: {cov_shape}")
    lines.append(f"Merged shape before subtype filter: {merged_before_shape}")
    lines.append(f"Merged shape after subtype filter: {merged_after_shape}")
    lines.append(f"Analysis shape after complete-case covariate filtering: {analysis_shape}")
    lines.append(f"Number of MSDL edges: {n_edges}\n")

    lines.append("[Subtype counts before complete-case covariate filtering]")
    lines.append(subtype_counts_before.to_string())
    lines.append("")

    lines.append("[Subtype counts after complete-case covariate filtering]")
    lines.append(subtype_counts_after.to_string())
    lines.append("")

    lines.append("[Base design matrix]")
    lines.append(f"Shape: {X_base_shape}")
    lines.append(f"Rank: {X_base_rank}")
    if dropped_rank_cols:
        lines.append("Dropped rank-deficient columns:")
        for c in dropped_rank_cols:
            lines.append(f"  - {c}")
    lines.append("")

    lines.append("[Global warnings]")
    if warnings_list:
        for w in warnings_list:
            lines.append(w)
    else:
        lines.append("No global feasibility warnings detected.")
    lines.append("")

    lines.append("[Edge diagnostics]")
    lines.append(f"Min non-missing subjects across edges: {edge_n_summary['min_n']}")
    lines.append(f"Median non-missing subjects across edges: {edge_n_summary['median_n']}")
    lines.append(f"Max non-missing subjects across edges: {edge_n_summary['max_n']}")
    lines.append(f"Near-zero variance edges before fitting: {near_zero_var_edges}")
    lines.append("")

    lines.append("[Fit summary]")
    lines.append(f"Edges successfully fit: {results_df.shape[0]}")
    lines.append(f"Edges skipped: {skipped_df.shape[0]}")
    lines.append("")

    if results_df.shape[0] > 0:
        lines.append("[Strict pattern counts]")
        lines.append(results_df["strict_pattern_label"].value_counts(dropna=False).to_string())
        lines.append("")
        lines.append("[Candidate pattern counts]")
        lines.append(results_df["candidate_pattern_label"].value_counts(dropna=False).to_string())
        lines.append("")

    if skipped_df.shape[0] > 0:
        lines.append("[Skipped edge reasons]")
        lines.append(skipped_df["reason"].value_counts(dropna=False).to_string())
        lines.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(lines))


# =========================================================
# ARGPARSE
# =========================================================
def parse_args():
    parser = argparse.ArgumentParser(description="Model 2: MSDL phenotype-only subtype model")
    parser.add_argument("--connectivity", type=str, default=str(DEFAULT_CONNECTIVITY_FILE))
    parser.add_argument("--covariate", type=str, default=str(DEFAULT_COVARIATE_FILE))
    parser.add_argument("--output_dir", type=str, default=str(DEFAULT_OUTPUT_DIR))
    return parser.parse_args()


# =========================================================
# MAIN
# =========================================================
def main():
    args = parse_args()

    connectivity_file = Path(args.connectivity)
    covariate_file = Path(args.covariate)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n======================")
    print("MODEL 2: MSDL PHENOTYPE-ONLY SUBTYPE MODEL")
    print("======================")
    print(f"Connectivity file : {connectivity_file}")
    print(f"Covariate file    : {covariate_file}")
    print(f"Output dir        : {output_dir}")

    if not connectivity_file.exists():
        raise FileNotFoundError(f"Connectivity file not found: {connectivity_file}")
    if not covariate_file.exists():
        raise FileNotFoundError(f"Covariate file not found: {covariate_file}")

    warnings_list = []

    # -------------------------
    # Load inputs
    # -------------------------
    conn = pd.read_csv(connectivity_file)
    cov = pd.read_csv(covariate_file)

    conn_shape = conn.shape
    cov_shape = cov.shape

    # -------------------------
    # Normalize subject IDs
    # -------------------------
    subject_col_conn = None
    for c in ["subject_id", "SubjectId", "subject"]:
        if c in conn.columns:
            subject_col_conn = c
            break
    if subject_col_conn is None:
        raise ValueError("Connectivity file must contain a subject_id column.")

    if "SubjectId" not in cov.columns:
        raise ValueError("Covariate file must contain SubjectId column.")

    conn = conn.copy()
    cov = cov.copy()

    conn["SubjectId"] = conn[subject_col_conn].map(normalize_subject_id)
    cov["SubjectId"] = cov["SubjectId"].map(normalize_subject_id)

    if "Subject_match_check" in cov.columns:
        before = cov.shape[0]
        cov = cov[cov["Subject_match_check"].fillna(True)].copy()
        after = cov.shape[0]
        print(f"[Covariate] Kept Subject_match_check=True rows: {after}/{before}")

    # -------------------------
    # Edge columns
    # -------------------------
    edge_cols = [c for c in conn.columns if c not in [subject_col_conn, "SubjectId"]]
    print("\n[Input shapes]")
    print(f"Connectivity raw shape : {conn.shape}")
    print(f"Covariate raw shape    : {cov.shape}")
    print(f"Connectivity edges     : {len(edge_cols)}")

    conn[edge_cols] = conn[edge_cols].apply(pd.to_numeric, errors="coerce")

    # -------------------------
    # Merge
    # -------------------------
    keep_cov_cols = [c for c in ["SubjectId", "Age", "Sex", "Scan_type", "Manufacturer", "Subtype"] if c in cov.columns]
    missing_cov = [c for c in ["Age", "Sex", "Scan_type", "Manufacturer", "Subtype"] if c not in keep_cov_cols]
    if missing_cov:
        raise ValueError(f"Covariate file missing required columns: {missing_cov}")

    cov_use = cov[keep_cov_cols].drop_duplicates(subset=["SubjectId"]).copy()
    merged_before = conn.merge(cov_use, on="SubjectId", how="inner")
    merged_before_shape = merged_before.shape

    merged = merged_before.copy()
    merged["Subtype_std"] = merged["Subtype"].map(normalize_subtype)
    merged = merged[merged["Subtype_std"].isin(["Control", "AsymAD", "TypicalAD"])].copy()
    merged_after_shape = merged.shape

    print("\n[Merged data]")
    print(f"Merged rows before subtype filter : {merged_before_shape[0]}")
    print(f"Merged rows after subtype filter  : {merged_after_shape[0]}")
    print(f"Merged shape                      : {merged_after_shape}")

    subtype_counts_before = merged["Subtype_std"].value_counts(dropna=False)
    print_counts_block("[Subtype counts before complete-case covariate filtering]", subtype_counts_before)

    # -------------------------
    # Build design matrix
    # -------------------------
    analysis_df, X_base = build_base_design(merged)
    X_base, dropped_rank_cols = reduce_to_full_rank(X_base)

    X_rank = np.linalg.matrix_rank(X_base.to_numpy(dtype=float))
    analysis_shape = analysis_df.shape

    print("\n[Complete-case covariate data]")
    print(f"Rows kept after covariate complete-case filtering: {analysis_df.shape[0]}")
    print(f"Base design matrix shape                         : {X_base.shape}")
    print(f"Base design matrix rank                          : {X_rank}")
    if dropped_rank_cols:
        print("[Design fix] Dropped rank-deficient columns:")
        for c in dropped_rank_cols:
            print(f"  - {c}")

    subtype_counts_after = analysis_df["Subtype_std"].value_counts(dropna=False)
    print_counts_block("[Subtype counts after complete-case covariate filtering]", subtype_counts_after)

    # -------------------------
    # Global diagnostics
    # -------------------------
    if dropped_rank_cols:
        warnings_list.append(
            "WARNING: Base design matrix was rank deficient; redundant columns were dropped: "
            + ", ".join(dropped_rank_cols)
        )

    min_group_n = subtype_counts_after.min() if len(subtype_counts_after) > 0 else 0
    if min_group_n < 10:
        warnings_list.append(
            f"WARNING: At least one subtype group has very small sample size (min n={min_group_n}). "
            "Interpret results as exploratory."
        )

    # -------------------------
    # Edge diagnostics before fit
    # -------------------------
    edge_nonmissing = analysis_df[edge_cols].notna().sum(axis=0)
    edge_variances = analysis_df[edge_cols].var(axis=0, ddof=1)
    near_zero_var_edges = int((edge_variances <= NEAR_ZERO_VAR).sum())

    edge_n_summary = {
        "min_n": int(edge_nonmissing.min()),
        "median_n": int(edge_nonmissing.median()),
        "max_n": int(edge_nonmissing.max()),
    }

    print("\n[Edge diagnostics before per-edge fitting]")
    print(f"Min non-missing subjects across edges : {edge_n_summary['min_n']}")
    print(f"Median non-missing subjects across edges : {edge_n_summary['median_n']}")
    print(f"Max non-missing subjects across edges : {edge_n_summary['max_n']}")
    print(f"Near-zero variance edges : {near_zero_var_edges}")

    print("\n[Warnings]")
    if warnings_list:
        for w in warnings_list:
            print(w)
    else:
        print("No global feasibility warnings detected.")

    merged_out = output_dir / "model2_input_merged_dataset.csv.gz"
    analysis_df.to_csv(merged_out, index=False, compression="gzip")

    # -------------------------
    # Fit all edges
    # -------------------------
    print("\n[Fitting edge models...]")
    results = []
    skipped = []
    subtype_series = analysis_df["Subtype_std"]

    for i, edge in enumerate(edge_cols, start=1):
        fit_out, skip_out = fit_edge_ols(X_base, analysis_df[edge], subtype_series)

        if fit_out is not None:
            fit_out["edge_id"] = edge
            results.append(fit_out)
        else:
            skip_out["edge_id"] = edge
            skipped.append(skip_out)

        if (i % 50 == 0) or (i == len(edge_cols)):
            print(f"  Processed edges: {i}/{len(edge_cols)}")

    results_df = pd.DataFrame(results)
    skipped_df = pd.DataFrame(skipped)

    if results_df.shape[0] == 0:
        diagnostic_out = output_dir / "model2_diagnostic_report.txt"
        skipped_out = output_dir / "model2_skipped_edges.csv"
        if skipped_df.shape[0] > 0:
            skipped_df.to_csv(skipped_out, index=False)

        write_diagnostic_report(
            out_path=diagnostic_out,
            connectivity_file=connectivity_file,
            covariate_file=covariate_file,
            output_dir=output_dir,
            conn_shape=conn_shape,
            cov_shape=cov_shape,
            merged_before_shape=merged_before_shape,
            merged_after_shape=merged_after_shape,
            analysis_shape=analysis_shape,
            n_edges=len(edge_cols),
            subtype_counts_before=subtype_counts_before,
            subtype_counts_after=subtype_counts_after,
            X_base_shape=X_base.shape,
            X_base_rank=X_rank,
            dropped_rank_cols=dropped_rank_cols,
            edge_n_summary=edge_n_summary,
            near_zero_var_edges=near_zero_var_edges,
            warnings_list=warnings_list,
            results_df=results_df,
            skipped_df=skipped_df,
        )
        raise RuntimeError("No valid edge models were fit. Check model2_diagnostic_report.txt and skipped edges file.")

    # -------------------------
    # Multiple testing + labels
    # -------------------------
    results_df = add_fdr_columns(results_df)
    results_df["strict_pattern_label"] = results_df.apply(classify_pattern_strict, axis=1)
    results_df["candidate_pattern_label"] = results_df.apply(classify_pattern_candidate, axis=1)

    results_df["is_strict_resilience"] = results_df["strict_pattern_label"].str.startswith("Resilience_", na=False)
    results_df["is_candidate_resilience"] = results_df["candidate_pattern_label"].str.startswith("CandidateResilience_", na=False)

    # Sort
    results_df = results_df.sort_values(
        ["fdr_subtype_omnibus", "p_subtype_omnibus", "p_asym_vs_typ", "p_typ_vs_control"],
        na_position="last"
    ).reset_index(drop=True)

    # -------------------------
    # Save outputs
    # -------------------------
    all_out = output_dir / "model2_all_edge_results.csv.gz"

    strict_res_out = output_dir / "model2_strict_resilience_edges.csv"
    strict_comp_out = output_dir / "model2_strict_compensation_edges.csv"
    strict_nonres_out = output_dir / "model2_strict_non_resilience_edges.csv"
    strict_mixed_out = output_dir / "model2_strict_mixed_unclear_edges.csv"

    cand_res_out = output_dir / "model2_candidate_resilience_edges.csv"
    cand_comp_out = output_dir / "model2_candidate_compensation_edges.csv"
    cand_nonres_out = output_dir / "model2_candidate_non_resilience_edges.csv"
    cand_mixed_out = output_dir / "model2_candidate_mixed_unclear_edges.csv"

    skipped_out = output_dir / "model2_skipped_edges.csv"
    top50_fdr_out = output_dir / "model2_top50_by_omnibus_fdr.csv"
    top50_raw_out = output_dir / "model2_top50_by_raw_typ_vs_control.csv"
    diagnostic_out = output_dir / "model2_diagnostic_report.txt"

    results_df.to_csv(all_out, index=False, compression="gzip")
    skipped_df.to_csv(skipped_out, index=False)

    # Strict subsets
    strict_resilience_df = results_df[results_df["strict_pattern_label"].str.startswith("Resilience_", na=False)].copy()
    strict_compensation_df = results_df[results_df["strict_pattern_label"].str.startswith("Compensation_", na=False)].copy()
    strict_nonres_df = results_df[results_df["strict_pattern_label"].str.startswith("NonResilience_", na=False)].copy()
    strict_mixed_df = results_df[
        ~(results_df["strict_pattern_label"].str.startswith("Resilience_", na=False)) &
        ~(results_df["strict_pattern_label"].str.startswith("Compensation_", na=False)) &
        ~(results_df["strict_pattern_label"].str.startswith("NonResilience_", na=False))
    ].copy()

    # Candidate subsets
    candidate_resilience_df = results_df[results_df["candidate_pattern_label"].str.startswith("CandidateResilience_", na=False)].copy()
    candidate_compensation_df = results_df[results_df["candidate_pattern_label"].str.startswith("CandidateCompensation_", na=False)].copy()
    candidate_nonres_df = results_df[results_df["candidate_pattern_label"].str.startswith("CandidateNonResilience_", na=False)].copy()
    candidate_mixed_df = results_df[
        ~(results_df["candidate_pattern_label"].str.startswith("CandidateResilience_", na=False)) &
        ~(results_df["candidate_pattern_label"].str.startswith("CandidateCompensation_", na=False)) &
        ~(results_df["candidate_pattern_label"].str.startswith("CandidateNonResilience_", na=False))
    ].copy()

    strict_resilience_df.to_csv(strict_res_out, index=False)
    strict_compensation_df.to_csv(strict_comp_out, index=False)
    strict_nonres_df.to_csv(strict_nonres_out, index=False)
    strict_mixed_df.to_csv(strict_mixed_out, index=False)

    candidate_resilience_df.to_csv(cand_res_out, index=False)
    candidate_compensation_df.to_csv(cand_comp_out, index=False)
    candidate_nonres_df.to_csv(cand_nonres_out, index=False)
    candidate_mixed_df.to_csv(cand_mixed_out, index=False)

    results_df.head(50).to_csv(top50_fdr_out, index=False)
    results_df.sort_values(["p_typ_vs_control", "p_asym_vs_typ", "p_asym_vs_control"]).head(50).to_csv(top50_raw_out, index=False)

    # -------------------------
    # Write diagnostic report
    # -------------------------
    write_diagnostic_report(
        out_path=diagnostic_out,
        connectivity_file=connectivity_file,
        covariate_file=covariate_file,
        output_dir=output_dir,
        conn_shape=conn_shape,
        cov_shape=cov_shape,
        merged_before_shape=merged_before_shape,
        merged_after_shape=merged_after_shape,
        analysis_shape=analysis_shape,
        n_edges=len(edge_cols),
        subtype_counts_before=subtype_counts_before,
        subtype_counts_after=subtype_counts_after,
        X_base_shape=X_base.shape,
        X_base_rank=X_rank,
        dropped_rank_cols=dropped_rank_cols,
        edge_n_summary=edge_n_summary,
        near_zero_var_edges=near_zero_var_edges,
        warnings_list=warnings_list,
        results_df=results_df,
        skipped_df=skipped_df,
    )

    # -------------------------
    # Final console summary
    # -------------------------
    print("\n======================")
    print("DONE")
    print("======================")
    print(f"All results                 : {all_out}")
    print(f"Strict resilience edges     : {strict_res_out} ({strict_resilience_df.shape[0]})")
    print(f"Candidate resilience edges  : {cand_res_out} ({candidate_resilience_df.shape[0]})")
    print(f"Strict compensation edges   : {strict_comp_out} ({strict_compensation_df.shape[0]})")
    print(f"Candidate compensation      : {cand_comp_out} ({candidate_compensation_df.shape[0]})")
    print(f"Strict non-resilience       : {strict_nonres_out} ({strict_nonres_df.shape[0]})")
    print(f"Candidate non-resilience    : {cand_nonres_out} ({candidate_nonres_df.shape[0]})")
    print(f"Strict mixed/unclear        : {strict_mixed_out} ({strict_mixed_df.shape[0]})")
    print(f"Candidate mixed/unclear     : {cand_mixed_out} ({candidate_mixed_df.shape[0]})")
    print(f"Skipped edges               : {skipped_out}")
    print(f"Diagnostic report           : {diagnostic_out}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()
