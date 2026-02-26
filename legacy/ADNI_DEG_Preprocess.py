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
import os

# --- 1. Define File Paths ---
EXPRESSION_FILE = 'ADNI_Gene_Expression_Profile.csv'
COVARIATES_FILE = 'covariates.csv'
OUTPUT_MATRIX_FILE = 'adni_expression_matrix_cleaned.csv'
OUTPUT_METADATA_FILE = 'adni_metadata_merged.csv'

print("--- Starting Pipeline ---")

# --- 2. Load the Raw Data ---
if not os.path.exists(EXPRESSION_FILE):
    print(f"Error: {EXPRESSION_FILE} not found.")
    exit()

print("Loading large expression file (this may take a moment)...")
# Read the file. Low_memory=False is needed because columns have mixed types (headers vs numbers)
data_full = pd.read_csv(EXPRESSION_FILE, header=None, low_memory=False)

# --- 3. Process Metadata (Top 8 Rows) ---
print("Processing metadata headers...")
# Extract rows 0-7, columns 3 onwards (where subject data starts)
metadata_raw = data_full.iloc[0:8, 3:].T.reset_index(drop=True)

# Rename columns manually based on ADNI structure
metadata_raw.columns = ["Phase", "Visit", "SubjectID_Original", "purity_260_280", 
                        "purity_260_230", "rin", "affy_qc", "globin_qc"]

# Create a 'SubjectID_Base' (e.g., 011_S_0010) for merging with covariates
# Assuming format is "PTID_VisitCode" or similar, usually the first 3 parts are the ID
metadata_raw['SubjectID_Base'] = metadata_raw['SubjectID_Original'].str.extract(r'(\d{3}_S_\d{4})')

# --- 4. Process Expression Data (Row 9 onwards) ---
print("Processing expression matrix...")
# Column 0 = ProbeSetID, Column 1 = LocusLink, Column 2 = GeneSymbol
# Data starts at Column 3
gene_info = data_full.iloc[8:, :3].copy()
gene_info.columns = ['ProbeSetID', 'LocusLink', 'GeneSymbol']

# Extract the numeric data
expression_values = data_full.iloc[8:, 3:].copy()

# Convert to numeric, forcing non-numeric to NaN (and fill with 0 or drop)
print("Converting expression values to numeric...")
expression_values = expression_values.apply(pd.to_numeric, errors='coerce')

# Add GeneSymbol to the data for collapsing
expression_values['GeneSymbol'] = gene_info['GeneSymbol'].values

# --- 5. Clean and Collapse Probes (The Fix) ---
print(f"Initial probe count: {len(expression_values)}")

# Remove rows with no Gene Symbol or empty strings
expression_values = expression_values.dropna(subset=['GeneSymbol'])
expression_values = expression_values[expression_values['GeneSymbol'] != '---']
expression_values = expression_values[expression_values['GeneSymbol'].str.strip() != '']

print(f"Count after removing missing symbols: {len(expression_values)}")

# Strategy: Max Mean (Keep the probe with highest average expression for each gene)
# 1. Calculate mean expression for each row (probe)
expression_values['Mean_Exp'] = expression_values.iloc[:, :-1].mean(axis=1)

# 2. Sort by Gene Symbol and Mean Expression (Ascending)
expression_values = expression_values.sort_values(by=['GeneSymbol', 'Mean_Exp'])

# 3. Drop duplicates keeping the LAST (which has the highest mean)
expression_values_collapsed = expression_values.drop_duplicates(subset=['GeneSymbol'], keep='last')

print(f"Final gene count (unique genes): {len(expression_values_collapsed)}")

# --- 6. Finalize Matrix Structure ---
# Set Index to GeneSymbol and remove the helper columns
final_matrix = expression_values_collapsed.set_index('GeneSymbol')
final_matrix = final_matrix.drop(columns=['Mean_Exp'])

# Assign Subject IDs as columns
final_matrix.columns = metadata_raw['SubjectID_Original']

# Transpose: Subjects as Rows, Genes as Columns (Standard for ML/Stats)
final_matrix_T = final_matrix.T

# --- 7. Merge with Clinical Covariates ---
print("Merging with clinical covariates...")
try:
    clinical_covariates = pd.read_csv(COVARIATES_FILE)
    
    # Ensure standard naming
    if 'subject_id' in clinical_covariates.columns:
        clinical_covariates.rename(columns={'subject_id': 'SubjectID_Base'}, inplace=True)
    
    # Merge
    # We use 'SubjectID_Base' to link the expression metadata with the clinical file
    full_metadata = pd.merge(metadata_raw, clinical_covariates, on='SubjectID_Base', how='inner')
    
    # Filter the expression matrix to only include subjects that exist in the clinical data
    # (The 'inner' merge might have dropped some subjects who didn't have clinical info)
    common_subjects = full_metadata['SubjectID_Original'].values
    final_matrix_T = final_matrix_T.loc[common_subjects]
    
    print(f"Merged successfully. Final Subject Count: {len(final_matrix_T)}")
    
except Exception as e:
    print(f"Warning: Issue merging covariates. Saving expression only. Error: {e}")
    full_metadata = metadata_raw

# --- 8. Save Outputs ---
print("Saving files...")
final_matrix_T.to_csv(OUTPUT_MATRIX_FILE)
full_metadata.to_csv(OUTPUT_METADATA_FILE, index=False)

# --- 9. Data Range Check ---
max_val = final_matrix_T.max().max()
print("\n--- Data Check ---")
print(f"Maximum expression value in dataset: {max_val:.2f}")
if max_val < 20:
    print("INSIGHT: Data appears to be Log2 transformed already (Standard for ADNI).")
else:
    print("INSIGHT: Data appears to be Raw Intensity. You MUST log-transform before DEG analysis.")

print("\nDone!")
