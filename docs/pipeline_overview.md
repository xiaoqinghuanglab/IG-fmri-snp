# Pipeline overview

This repository organizes code for a subtype-aware imaging-genetics workflow in AD.

**Main stages**
1. **GWAS (PLINK):** run logistic regression per subtype contrast and keep the top-ranked SNPs per contrast.
   Output: `snp_metadata.csv` + `genotype_matrix.csv` (used as SNP features).
2. **rs-fMRI preprocessing (FSL):** motion correction, slice timing, MNI alignment, smoothing.
   Output: subject-level preprocessed 4D NIfTIs.
3. **Connectivity features (Nilearn + scikit-learn):** extract MSDL ROI time series and compute sparse connectivity via GraphicalLassoCV.
   Output: `Connectivity_matrix_all_subjects_region_pairs.csv`
4. **Imaging-genetics association (statsmodels):** OLS with SNP×Subtype interactions for each SNP-edge pair + BH-FDR.
   Output: `all_regression_result*.csv`, `significant_result.csv`
5. **Cross-cohort validation (blood expression + MRI):** scripts in `scripts/transcriptomics/` and `scripts/validation/`.

Each stage can be run independently and all scripts accept **input/output paths as CLI flags** (no code edits).
