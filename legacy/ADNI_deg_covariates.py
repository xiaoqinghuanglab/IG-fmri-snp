"""Legacy (original) analysis script.

This file is kept for provenance and for readers who want to see the original,
notebook-style workflow. For a fully parameterized version (recommended for
reproducible runs), see the corresponding script under:

  scripts/transcriptomics/

Notes
-----
- Paths in this legacy script are hard-coded (as in the original analysis).
- New users should prefer the parameterized scripts + README instructions.
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import time
import sys

# --- Configuration ---
MATRIX_FILE = 'adni_expression_matrix_cleaned.csv'
METADATA_FILE = 'adni_metadata_merged.csv'

# Define the exact group names as found in your "Subtypes" column
# Format: (Test_Group, Reference_Group)
COMPARISONS = [
    ('Asym AD', 'Typical AD'),
    ('Asym AD', 'Control'),
    ('Typical AD', 'Control'),
    ('Asym AD', 'Low-NFT AD'),
    ('Low-NFT AD', 'Typical AD'),
    ('Low-NFT AD', 'Control')
]

print("--- Loading Data ---")
# Load data (Subjects are Rows, Genes are Columns)
expression_data = pd.read_csv(MATRIX_FILE, index_col=0) 
metadata = pd.read_csv(METADATA_FILE)

print(f"Expression Matrix Shape: {expression_data.shape} (Subjects x Genes)")

# --- Preprocessing Covariates ---
print("Cleaning Covariates...")

# 1. Encode Gender (M/F -> 0/1)
# Check unique values first to be safe
if metadata['Gender'].dtype == object:
    metadata['Gender_Code'] = metadata['Gender'].astype('category').cat.codes
else:
    metadata['Gender_Code'] = metadata['Gender']

# 2. Ensure APOE4 is numeric
metadata['APOE4'] = pd.to_numeric(metadata['APOE4'], errors='coerce')

# 3. Handle Missing Covariates
required_cols = ['SubjectID_Original', 'Subtypes', 'Age', 'Gender_Code', 'APOE4']
metadata_clean = metadata.dropna(subset=required_cols).copy()
print(f"Subjects with full covariate data: {len(metadata_clean)}")

# --- Define the Regression Function ---
def run_deg_analysis(test_group, ref_group, meta_df, expr_df):
    print(f"\n>>> Running Comparison: {test_group} vs {ref_group}")
    
    # 1. Filter Metadata for this pair
    sub_meta = meta_df[meta_df['Subtypes'].isin([test_group, ref_group])].copy()
    
    # 2. Encode Group (Reference = 0, Test = 1)
    sub_meta['Group_Code'] = np.where(sub_meta['Subtypes'] == test_group, 1, 0)
    
    # 3. Sync Subjects between Metadata and Expression Data
    # Get subjects that exist in BOTH the metadata subset AND the expression matrix index
    target_subs = sub_meta['SubjectID_Original'].values
    valid_subjects = [s for s in target_subs if s in expr_df.index]
    
    # Check sample size
    n_test = sub_meta[sub_meta['Subtypes'] == test_group]['SubjectID_Original'].isin(valid_subjects).sum()
    n_ref = sub_meta[sub_meta['Subtypes'] == ref_group]['SubjectID_Original'].isin(valid_subjects).sum()
    
    print(f"   Subjects found -> {test_group}: {n_test}, {ref_group}: {n_ref}")
    
    if n_test < 3 or n_ref < 3:
        print("   SKIP: Not enough subjects in one or both groups.")
        return

    # 4. Prepare Matrices for Regression
    # Align Metadata and Expression by Subject ID
    current_meta = sub_meta.set_index('SubjectID_Original').loc[valid_subjects]
    current_expr = expr_df.loc[valid_subjects] # Rows=Subjects, Cols=Genes
    
    # Design Matrix (X): Intercept + Group + Age + Gender + APOE4
    X = current_meta[['Group_Code', 'Age', 'Gender_Code', 'APOE4']]
    X = sm.add_constant(X)
    
    # Data Matrix (Y)
    Y_matrix = current_expr.values
    gene_names = current_expr.columns
    total_genes = len(gene_names)
    
    print(f"   Fitting models for {total_genes} genes...")
    
    # 5. Loop through genes (Vectorized preparation)
    results_list = []
    start_time = time.time()
    
    # Pre-calculate X inverse stuff to speed up OLS? 
    # For simplicity and reliability, we stick to the loop but use simple OLS.
    X_val = X.values
    
    for i in range(total_genes):
        y = Y_matrix[:, i]
        
        # OLS Regression
        try:
            model = sm.OLS(y, X_val)
            result = model.fit()
            
            # Index 1 is 'Group_Code' (our variable of interest)
            results_list.append({
                'GeneSymbol': gene_names[i],
                'Log2FC': result.params[1],
                'P_Value': result.pvalues[1],
                'Age_P': result.pvalues[2] if len(result.pvalues) > 2 else np.nan
            })
        except:
            continue

        # Progress bar every 2000 genes
        if i % 2000 == 0 and i > 0:
            sys.stdout.write(f"\r   Progress: {i}/{total_genes} genes")
            sys.stdout.flush()
            
    print(f"\r   Comparison finished in {time.time() - start_time:.1f} seconds.      ")
    
    # 6. FDR Correction & Saving
    res_df = pd.DataFrame(results_list)
    
    if res_df.empty:
        print("   Error: No results generated.")
        return

    # Remove NaNs
    res_df = res_df.dropna(subset=['P_Value'])
    
    # FDR
    reject, pvals_corrected, _, _ = multipletests(res_df['P_Value'], alpha=0.05, method='fdr_bh')
    res_df['FDR'] = pvals_corrected
    
    # Save
    res_df = res_df.sort_values(by='P_Value')
    filename = f"DEG_{test_group.replace(' ', '')}_vs_{ref_group.replace(' ', '')}_Covariates.csv"
    res_df.to_csv(filename, index=False)
    
    sig_count = len(res_df[res_df['FDR'] < 0.05])
    print(f"   Saved to {filename}")
    print(f"   Significant Genes (FDR < 0.05): {sig_count}")

# --- Run All Comparisons ---
for test, ref in COMPARISONS:
    run_deg_analysis(test, ref, metadata_clean, expression_data)

print("\n--- All Analyses Complete ---")
