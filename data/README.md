# Bundled data inputs

This repository keeps only the minimal CSV inputs needed to reproduce the paper's downstream analysis after cloning.

## Files intentionally kept in `data/`

- `covariate_file.csv`
- `msdl_all_subjects_connectivity_edges.csv`
- `genotype_matrix.csv`
- `snp_metadata.csv`

These are enough to:

1. generate candidate resilience edges with
   `bash run.sh python scripts/fmriphenotype/get_resilience_edge.py`
2. run the pairwise SNP-connectivity analysis with
   `bash run.sh imggenetics`

## Files intentionally not kept in `data/`

- `model2_candidate_resilience_edges.csv`
  This is a generated output from the phenotype-only stage and is written to `outputs/fmriphenotype/`.
- per-contrast intermediate genotype/SNP files
  These are upstream GWAS intermediates and are replaced here by the final `genotype_matrix.csv` and `snp_metadata.csv`.
- alternate or outdated connectivity / covariate files
  The cleaned paper workflow is MSDL-only and uses the current `covariate_file.csv`.

## When you need upstream regeneration

If you want to rebuild the bundled inputs from raw data instead of using the included CSVs:

- connectivity extraction instructions are in [scripts/fmri/README.md](/Users/shihab/Codex/IG-fmri-snp/scripts/fmri/README.md:1)
- GWAS preparation instructions are in [scripts/gwas/README.md](/Users/shihab/Codex/IG-fmri-snp/scripts/gwas/README.md:1)

Raw ADNI / ANMerge data, raw PLINK files, and raw fMRIPrep derivatives are not distributed in this repo.
