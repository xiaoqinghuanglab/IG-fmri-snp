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
├── fmriprep_derivatives/
│   └── sub-XXXX/func/...          # input to MSDL connectivity extraction
├── imggenetics_inputs/
│   ├── genotype_matrix.csv
│   ├── snp_metadata.csv
│   ├── covariate_file.csv
│   ├── msdl_all_subjects_connectivity_edges.csv
│   └── model2_candidate_resilience_edges.csv
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
If you only want to run the pairwise SNP-connectivity analysis, you need **five CSV files**:

- `genotype_matrix.csv`
- `snp_metadata.csv`
- `covariate_file.csv`
- `msdl_all_subjects_connectivity_edges.csv`
- `model2_candidate_resilience_edges.csv`

You can pass these as explicit paths using:
`--genetic-file`, `--connectivity-file`, `--covariates-file`, `--snp-metadata-file`, and `--candidate-edges-file`.
