# Pipeline Overview

This repository implements a subtype-aware imaging-genetics workflow for Alzheimer's disease.

## Main Stages

1. **GWAS feature prioritization (PLINK)**
   Run pairwise subtype GWAS analyses and keep the top-ranked SNPs per contrast.
   Output: `genotype_matrix.csv` and `snp_metadata.csv`

2. **rs-fMRI preprocessing (fMRIPrep)**
   Generate MNI-space preprocessed resting-state BOLD files plus confounds tables.
   Output: fMRIPrep derivatives under a BIDS derivatives directory

3. **Connectivity extraction (MSDL + GraphicalLassoCV)**
   Extract MSDL ROI time series, regress confounds, apply filtering and z-scoring, and compute signed partial-correlation edges.
   Output: `msdl_all_subjects_connectivity_edges.csv`

4. **Phenotype-only subtype modeling**
   Identify candidate resilience-related connectivity edges before SNP-level interpretation.
   Output: `model2_candidate_resilience_edges.csv` plus related diagnostic tables

5. **Pairwise SNP-connectivity interaction analysis**
   For `AsymAD_vs_TypAD`, `TypAD_vs_Control`, and `AsymAD_vs_Control`, fit:
   `edge ~ SNP_c + Group + SNP_c:Group + age + sex + scan_type + manufacturer`
   Output: per-comparison full tables, significant hits, and lead-SNP-per-edge summaries

6. **Biological interpretation**
   Export VEP rsID lists, BrainNetViewer edge files, and optional transcriptomic / MRI follow-up analyses.

All runnable stages accept explicit input and output paths so the workflow can be reproduced without editing source files.
