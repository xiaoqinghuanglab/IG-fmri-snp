# Data layout (local only)

This repo does **not** include ADNI/ANMerge protected datasets.

You can keep your data anywhere on your machine and point scripts to it using `--input-*` and `--output-*` flags.

## Recommended local folder layout (example)

```text
my_data_root/
├── genotypes/
│   ├── QC_GENOTYPES_PREFIX.bed
│   ├── QC_GENOTYPES_PREFIX.bim
│   └── QC_GENOTYPES_PREFIX.fam
├── fmri/
│   ├── raw_nifti/                 # input to FSL preprocessing
│   └── processed_nifti/           # output from FSL preprocessing
├── imggenetics_inputs/
│   ├── genotype_matrix.csv
│   ├── snp_metadata.csv
│   ├── subject_covariates.csv
│   └── Connectivity_matrix_all_subjects_region_pairs.csv
├── transcriptomics/
│   ├── adni/
│   │   ├── ADNI_Gene_Expression_Profile.csv
│   │   └── covariates.csv
│   └── anmerge/
│       ├── ANMerge_Unified_Ready_for_DEA.csv
│       ├── ANMerge_MRI_FS6.0_under_90.csv
│       ├── ANMerge_NewModel_NewLabel_Clinical_Info.csv
│       └── ANMerge_blood_rna_gene_expr_removedbatch_02022024.csv
└── outputs/                        # recommended to keep outputs outside the repo
```

## Minimal inputs for `imggenetics`
If you only want to run the SNP×Subtype OLS, you need **four CSV files**:

- `genotype_matrix.csv` (must include `Subject_ID`)
- `Connectivity_matrix_all_subjects_region_pairs.csv` (must include `Subject_ID`)
- `subject_covariates.csv` (must include `Subject_ID`, `Subtype`, `age`, `sex`, `scanner`, `scan_type`)
- `snp_metadata.csv` (must include `SNP`, `Case`, `Control`)

You can pass these as explicit paths using:
`--genetic-file`, `--connectivity-file`, `--covariates-file`, `--snp-metadata-file`.
