# GWAS preparation workflow

Use the scripts in this folder when you want to regenerate the GWAS-derived SNP inputs instead of using the bundled `data/genotype_matrix.csv` and `data/snp_metadata.csv`.

## Inputs

- QC-passed PLINK prefix: `<PREFIX>.bed/.bim/.fam`
- a clinical / covariate table with subject IDs, subtype labels, age, sex, APOE4, and ancestry PCs

Only the three manuscript groups are used in the cleaned workflow:

- `Control`
- `AsymAD`
- `TypAD`

## Typical sequence

1. Build PLINK phenotype and covariate files:

```bash
bash run.sh python scripts/gwas/prepare_pheno_covar.py \
  --clinical /path/to/clinical_covariates.csv \
  --out-dir outputs/gwas/pheno_covar
```

2. Run contrast-wise PLINK logistic regression:

```bash
bash run.sh bash scripts/gwas/run_gwas_plink.sh \
  --bfile /path/to/QC_GENOTYPES_PREFIX \
  --pheno-dir outputs/gwas/pheno_covar \
  --out-dir outputs/gwas/assoc
```

3. Extract top SNPs per contrast:

```bash
bash run.sh python scripts/gwas/extract_top_snps.py \
  --assoc outputs/gwas/assoc/gwas_AsymAD_vs_Control.assoc.logistic \
  --top-n 1000 \
  --out outputs/gwas/top1000_AsymAD_vs_Control.tsv
```

Repeat step 3 for:

- `AsymAD_vs_TypAD`
- `TypAD_vs_Control`

4. Build the final downstream SNP inputs:

```bash
bash run.sh python scripts/gwas/make_genotype_matrix.py \
  --bfile /path/to/QC_GENOTYPES_PREFIX \
  --top-list outputs/gwas/top1000_AsymAD_vs_Control.tsv AsymAD Control \
  --top-list outputs/gwas/top1000_AsymAD_vs_TypAD.tsv AsymAD TypAD \
  --top-list outputs/gwas/top1000_TypAD_vs_Control.tsv TypAD Control \
  --out-dir outputs/gwas/final
```

## Final outputs used by `imggenetics`

- `genotype_matrix.csv`
- `snp_metadata.csv`

You can either:

- pass them directly to `imggenetics` with explicit flags, or
- copy them into the repo `data/` folder if you are replacing the bundled reproduction inputs.
