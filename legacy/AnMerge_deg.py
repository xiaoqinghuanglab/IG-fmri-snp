import pandas as pd
import numpy as np
from scipy import stats
import os
import sys

# --- CONFIGURATION ---
INPUT_FILE = 'ANMerge_Unified_Ready_for_DEA.csv'
OUTPUT_DIR = 'DEG_Results_CovariateAdjusted'  # NEW folder
DX_COL = 'subtypes'

# Covariates to adjust for (must exist in the CSV)
COVARIATES = ['Age', 'sex']

COMPARISONS = [
    ('Asym AD', 'Typical AD'),
    ('Asym AD', 'Control'),
    ('Asym AD', 'Low-NFT AD'),
    ('Low-NFT AD', 'Typical AD'),
    ('Low-NFT AD', 'Control'),
    ('Typical AD', 'Control')
]

# --- 1. SETUP ---
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"--- Loading Data: {INPUT_FILE} ---")
try:
    df = pd.read_csv(INPUT_FILE)
    print(f"Data Loaded. Shape: {df.shape}")
except FileNotFoundError:
    print(f"[ERROR] Could not find {INPUT_FILE}.")
    sys.exit(1)

# --- 2. BASIC VALIDATION ---
missing_cols = [c for c in [DX_COL] + COVARIATES if c not in df.columns]
if missing_cols:
    print(f"[ERROR] Missing required columns in input file: {missing_cols}")
    sys.exit(1)

# --- 3. IDENTIFY GENES ---
# Exclude known clinical columns + covariates from gene list
clinical_cols = ['Subject_Id', DX_COL] + COVARIATES
gene_cols = [c for c in df.columns if c not in clinical_cols]
print(f"Identified {len(gene_cols)} genes for analysis.")

# --- 4. HELPERS ---
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR. Returns q-values aligned with pvals order."""
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

def build_design_matrix(sub: pd.DataFrame, group1: str, covariates: list[str]):
    """
    Builds X with columns:
      Intercept, group_indicator (1 if group1 else 0), covariates (dummy-coded if categorical)
    Returns: X (np.ndarray), colnames (list[str]), group_indicator (np.ndarray), row_mask (np.ndarray)
    """
    # Group indicator: 1 for group1, 0 for group2
    group_ind = (sub[DX_COL] == group1).astype(int)

    cov_df = sub[covariates].copy()

    processed_cols = []
    for col in cov_df.columns:
        s = cov_df[col]

        # If numeric-ish, coerce to numeric
        if pd.api.types.is_numeric_dtype(s):
            cov_df[col] = pd.to_numeric(s, errors='coerce')
            processed_cols.append(col)
        else:
            # Treat as categorical -> one-hot, drop_first to avoid collinearity
            cov_df[col] = s.astype('category')

    # One-hot encode categorical columns
    cov_df = pd.get_dummies(cov_df, drop_first=True)

    X_df = pd.concat([group_ind.rename('__group__'), cov_df], axis=1)

    # Add intercept
    X = np.column_stack([np.ones(len(X_df), dtype=float), X_df.to_numpy(dtype=float)])
    colnames = ['Intercept'] + list(X_df.columns)

    # Row mask for covariate completeness
    row_mask = np.isfinite(X).all(axis=1)

    return X, colnames, group_ind.to_numpy(dtype=int), row_mask

# --- 5. ANALYSIS FUNCTION (Covariate-adjusted OLS) ---
def run_covariate_deg(data: pd.DataFrame, group1: str, group2: str, genes: list[str], covariates: list[str]):
    print(f"\nProcessing (Covariate-adjusted OLS): {group1} vs {group2} ...")

    sub = data[data[DX_COL].isin([group1, group2])].copy()
    n_total = sub.shape[0]
    n1 = (sub[DX_COL] == group1).sum()
    n2 = (sub[DX_COL] == group2).sum()
    print(f" -> Samples: {group1}={n1}, {group2}={n2} (Total={n_total})")

    if n1 < 3 or n2 < 3:
        print(" -> [WARNING] Not enough samples for covariate-adjusted regression. Skipping.")
        return None

    # Build design matrix once per comparison
    X, xcols, group_ind, base_mask = build_design_matrix(sub, group1, covariates)
    p = X.shape[1]  # number of parameters

    results = []

    for i, gene in enumerate(genes):
        y = pd.to_numeric(sub[gene], errors='coerce').to_numpy(dtype=float)

        mask = base_mask & np.isfinite(y)
        n = int(mask.sum())

        # Need enough dof
        if n <= p + 1:
            continue

        Xm = X[mask, :]
        ym = y[mask]

        # Check rank (avoid singular design)
        if np.linalg.matrix_rank(Xm) < p:
            continue

        # OLS via least squares
        beta, _, _, _ = np.linalg.lstsq(Xm, ym, rcond=None)
        resid = ym - (Xm @ beta)

        df_resid = n - p
        if df_resid <= 0:
            continue

        s2 = (resid @ resid) / df_resid
        XtX_inv = np.linalg.inv(Xm.T @ Xm)
        se = np.sqrt(np.diag(s2 * XtX_inv))

        # Group effect is column "__group__" which is xcols[1]
        group_idx = 1
        if se[group_idx] == 0 or not np.isfinite(se[group_idx]):
            continue

        t_stat = beta[group_idx] / se[group_idx]
        p_val = 2.0 * stats.t.sf(np.abs(t_stat), df=df_resid)

        # Unadjusted means for reference (on the same masked rows)
        gmask = mask  # valid rows for this gene
        y1 = y[gmask & (group_ind == 1)]
        y2 = y[gmask & (group_ind == 0)]
        mean1 = float(np.nanmean(y1)) if y1.size else np.nan
        mean2 = float(np.nanmean(y2)) if y2.size else np.nan

        results.append({
            'Gene': gene,
            'Log2FC_adj': float(beta[group_idx]),   # adjusted group difference (group1 - group2)
            'P_Value': float(p_val),
            'T_Stat': float(t_stat),
            'N_Used': int(n),
            'Mean_' + group1.replace(' ', ''): mean1,
            'Mean_' + group2.replace(' ', ''): mean2
        })

        if (i + 1) % 1000 == 0:
            print(f"    ... processed {i + 1}/{len(genes)} genes")

    if not results:
        print(" -> [WARNING] No valid results generated.")
        return None

    res_df = pd.DataFrame(results)

    # Drop invalid p-values
    res_df = res_df.dropna(subset=['P_Value'])
    res_df = res_df[res_df['P_Value'].between(0, 1, inclusive="both")]

    if res_df.empty:
        print(" -> [WARNING] All results invalid after filtering.")
        return None

    # FDR correction
    res_df['FDR'] = bh_fdr(res_df['P_Value'].to_numpy())
    res_df.sort_values('P_Value', inplace=True)

    return res_df

# --- 6. EXECUTE LOOPS ---
for g1, g2 in COMPARISONS:
    deg_res = run_covariate_deg(df, g1, g2, gene_cols, COVARIATES)

    if deg_res is not None:
        clean_g1 = g1.replace(' ', '').replace('-', '')
        clean_g2 = g2.replace(' ', '').replace('-', '')
        filename = f"{clean_g1}_vs_{clean_g2}_CovAdj_OLS_DEG.csv"
        filepath = os.path.join(OUTPUT_DIR, filename)

        deg_res.to_csv(filepath, index=False)

        sig_count = int((deg_res['FDR'] < 0.05).sum())
        print(f" -> SAVED: {filename}")
        print(f" -> Significant Genes (FDR < 0.05): {sig_count}")

print(f"\nAll analyses finished. Check '{OUTPUT_DIR}' folder.")

