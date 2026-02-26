# GWAS (PLINK) workflow for subtype contrasts

This repository treats GWAS as a **feature prioritization** step: for each subtype contrast, variants are ranked by association p-value and the top set is carried forward into SNP×connectivity modeling.

The original, step-by-step commands used in this project are preserved as a PDF in:
- `docs/GWAS_Analysis_Documentation.pdf`

This document provides a repo-style overview of the same workflow.

---

## Inputs

### Genotypes
- QC-ready PLINK files: `<PREFIX>.bed`, `<PREFIX>.bim`, `<PREFIX>.fam`

### Clinical / covariates table
A per-subject table containing:
- subject identifier matching the PLINK `.fam`
- subtype labels (e.g., Control / Asym AD / Typical AD)
- covariates (age, sex, APOE ε4 status, ancestry PCs)

---

## Outputs (used downstream)

### 1) `snp_metadata.csv`
A long table that records (at minimum):
- the SNP rsID
- which contrast it came from (`Case`, `Control`)
- p-values and/or summary stats (optional but useful)

### 2) `genotype_matrix.csv`
A subject-by-SNP matrix used by `imggenetics`:
- rows: subjects (`Subject_ID`)
- columns: SNPs encoded as allele-count variables (0/1/2)
- values: additive genotype coding

---

## Suggested workflow

### A) QC (example checklist)
The original workflow applied standard filters such as:
- per-subject and per-variant missingness thresholds
- minor allele frequency thresholding
- Hardy–Weinberg filtering in controls
- removal of strand-ambiguous A/T and C/G variants
- ancestry PCA to derive PCs for covariate adjustment

Exact thresholds and commands depend on your cohort, platform, and project policy—use the included PDF as a concrete starting point.

### B) Association testing (per contrast)
Use PLINK logistic regression with additive coding and covariate adjustment (age/sex/APOE4/PC1–PC5). Run separately for each contrast (e.g., AsymAD vs Control, AsymAD vs TypAD, TypAD vs Control).

### C) Post-GWAS feature selection
For each contrast:
1. sort by p-value,
2. keep the top N variants (e.g., 1,000),
3. take the union across contrasts,
4. export the final genotype feature matrix and a metadata table.

---

## Helper scripts in this repo

This repo includes small helpers that make the above steps easier to reproduce:

- `scripts/gwas/prepare_pheno_covar.py`  
  Build per-contrast phenotype + covariate files for PLINK.

- `scripts/gwas/run_gwas_plink.sh`  
  Template runner for contrast-wise PLINK logistic regression.

- `scripts/gwas/extract_top_snps.py`  
  Extract top-N hits per contrast from `*.assoc.logistic` outputs.

- `scripts/gwas/make_genotype_matrix.py`  
  Build `genotype_matrix.csv` + `snp_metadata.csv` from top-hit lists.

For end-to-end example commands, see the GWAS section in the root `README.md`.
