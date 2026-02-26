# IG-fmri-snp — subtype-aware imaging genetics for rs-fMRI connectivity

This repository collects a complete **imaging-genetics workflow** that connects:

- **Genetics:** GWAS-based SNP feature prioritization per subtype contrast (PLINK)
- **Imaging:** rs-fMRI preprocessing (FSL) and connectivity feature extraction (MSDL 39 ROI → 741 edges)
- **Association testing:** SNP × subtype interaction models for SNP–edge pairs (Python / statsmodels)
- **Downstream interpretation:** VEP rsID lists, BrainNetViewer edge exports, enrichment/PPI helpers
- **External validation (optional):** blood transcriptomics (ADNI + ANMerge) and ANMerge expression × MRI models

All runnable steps are exposed as **command-line scripts** with **explicit input/output path flags** so users can run the same code without editing source files.

> **Protected data:** ADNI/ANMerge datasets are not distributed with this repo. See `data/README.md` for an example local layout and required filenames/columns.

---

## Repository layout

```text
.
├── src/imggenetics/                 # installable package: SNP×Subtype OLS for connectivity phenotypes
├── scripts/
│   ├── fmri/                        # FSL preprocessing + MSDL connectivity derivation
│   ├── gwas/                        # PLINK helpers + top-SNP extraction + genotype matrix builder
│   ├── imggenetics/                 # post-processing utilities (VEP lists, BrainNetViewer exports)
│   ├── transcriptomics/             # ADNI + ANMerge blood DEG analysis
│   └── validation/                  # ANMerge expression × MRI interaction models
├── notebooks/                       # original notebook version (reference)
├── docs/                            # additional workflow notes + GWAS documentation PDF
├── data/README.md                   # data layout (local only)
├── requirements.txt                 # Python dependencies
└── run.sh                           # convenience runner (creates .venv, installs deps, runs commands)
```

---

## Installation

```bash
git clone https://github.com/xiaoqinghuanglab/IG-fmri-snp.git
cd IG-fmri-snp

# recommended: use the wrapper (creates .venv + installs deps)
bash run.sh python -c "import imggenetics; print('ok')"
```

You can also install manually:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

---

## Quick start: run the imaging-genetics OLS (SNP × subtype)

### What this step does
Given:
- a **genotype feature matrix** (subjects × SNPs),
- an **edge feature matrix** (subjects × 741 connectivity edges),
- **subject covariates** (age/sex/scanner/scan_type + subtype label),
- and a **SNP metadata table** that groups SNPs by their originating GWAS contrast,

this stage fits an OLS model per SNP–edge pair with a genotype-by-subtype interaction and writes:
- full model results
- BH-FDR corrected p-values
- a filtered “significant” table

### Required inputs
You can place these 4 CSVs anywhere and pass explicit paths:

| File | Required columns (minimum) |
|---|---|
| `genotype_matrix.csv` | `Subject_ID` + SNP columns (e.g., `rs12345_A`) |
| `Connectivity_matrix_all_subjects_region_pairs.csv` | `Subject_ID` + connectivity columns |
| `subject_covariates.csv` | `Subject_ID`, `Subtype`, `age`, `sex`, `scanner`, `scan_type` |
| `snp_metadata.csv` | `SNP` (rsID), `Case`, `Control` |

### Run
```bash
bash run.sh imggenetics   --output-dir outputs/imggenetics   --genetic-file /path/to/genotype_matrix.csv   --connectivity-file /path/to/Connectivity_matrix_all_subjects_region_pairs.csv   --covariates-file /path/to/subject_covariates.csv   --snp-metadata-file /path/to/snp_metadata.csv   --subtypes Control TypAD AsymAD   --reference-subtype Control
```

### Outputs
Written under `--output-dir`:
- `all_regression_result.csv`
- `all_regression_result_FDR.csv`
- `significant_result.csv`

---

## Stage 1 — rs-fMRI preprocessing (FSL)

### What this step does
A minimal rs-fMRI preprocessing routine:
- motion correction
- slice timing correction
- mean-image registration to MNI space
- application of the transform to the full 4D series
- smoothing

### Inputs
- A folder of **4D NIfTI** files (`*.nii.gz`) per subject/session.

### Outputs
- One preprocessed 4D NIfTI per input file, default naming: `<base>_smooth.nii.gz`.

### Run
```bash
bash run.sh bash scripts/fmri/preprocess_fmri_fsl.sh   --input-dir /path/to/raw_nifti   --output-dir /path/to/processed_nifti   --tr 3   --smooth-sigma 2
```

---

## Stage 2 — connectivity features (MSDL + GraphicalLassoCV)

### What this step does
For each subject’s preprocessed 4D NIfTI:
- extracts ROI time series using the **MSDL atlas (39 ROIs)**,
- optionally aligns time-series length (truncate/pad),
- fits **GraphicalLassoCV** and exports the (optionally absolute) precision edges.

This produces **741** edge features per subject (39×38/2).

### Inputs
- A folder of preprocessed NIfTIs (default pattern: `*_smooth.nii.gz`).
- Filenames should contain the subject identifier (default regex: `\d{3}_S_\d{4}`), or provide `--id-regex`.

### Outputs
- A subject-by-edge CSV (default column names: `Connectivity_<ROI1>_<ROI2>`).

### Run
```bash
bash run.sh python scripts/fmri/compute_connectivity_msdl.py   --input-dir /path/to/processed_nifti   --pattern "*_smooth.nii.gz"   --output-csv outputs/connectivity/Connectivity_matrix_all_subjects_region_pairs.csv   --tr 3   --fixed-length median   --abs-edges
```

---

## Stage 3 — GWAS for SNP feature prioritization (PLINK)

### What this stage is for
In this project, GWAS is used as a **feature selection step**: per contrast, SNPs are ranked by association p-value and the top set is carried forward into the SNP×connectivity modeling.

A full, command-oriented GWAS walkthrough is included in:
- `docs/GWAS_Analysis_Documentation.pdf` (original notes)
- `docs/gwas_workflow.md` (repo-style paraphrase)

### Typical inputs
- QC-passed PLINK files: `<PREFIX>.bed/.bim/.fam`
- a clinical/covariate table with: subject IDs, subtype labels, age/sex/APOE4, and ancestry PCs (PC1–PC5)

### Typical outputs
- contrast-specific PLINK association results (`*.assoc.logistic`)
- `top1000_<contrast>.tsv` files
- `snp_metadata.csv` (contrast labels + top SNPs)
- `genotype_matrix.csv` (subjects × unique top SNPs; 0/1/2 allele counts)

### Minimal run (helpers in `scripts/gwas/`)
```bash
# 1) build per-contrast phenotype + covariate files
bash run.sh python scripts/gwas/prepare_pheno_covar.py   --clinical-csv /path/to/clinical_covariates.csv   --out-dir outputs/gwas/pheno_covar

# 2) run PLINK logistic regression per contrast (template)
bash run.sh bash scripts/gwas/run_gwas_plink.sh   --bfile /path/to/QC_GENOTYPES_PREFIX   --pheno-dir outputs/gwas/pheno_covar   --out-dir outputs/gwas/assoc

# 3) extract top SNPs per contrast
bash run.sh python scripts/gwas/extract_top_snps.py   --assoc-dir outputs/gwas/assoc   --out-dir outputs/gwas   --top-n 1000

# 4) build genotype_matrix.csv + snp_metadata.csv for imggenetics
bash run.sh python scripts/gwas/make_genotype_matrix.py   --bfile /path/to/QC_GENOTYPES_PREFIX   --top-list outputs/gwas/top1000_AsymAD_vs_Control.tsv AsymAD Control   --top-list outputs/gwas/top1000_AsymAD_vs_TypAD.tsv AsymAD TypAD   --top-list outputs/gwas/top1000_TypAD_vs_Control.tsv TypAD Control   --out-dir outputs/gwas/final
```

---

## Stage 4 — post-processing (VEP lists + BrainNetViewer edge files)

These are optional utilities that turn `imggenetics` outputs into downstream-friendly files.

### VEP rsID lists + per-contrast result CSVs
**Input:** `significant_result.csv`  
**Output:** `VEP_Input_Files/*.txt` + `GWAS_Results_by_Comparison/*_results.csv`

```bash
bash run.sh python scripts/imggenetics/export_vep_inputs.py   --significant-csv outputs/imggenetics/significant_result.csv   --out-dir outputs/imggenetics_post
```

### BrainNetViewer edge exports
**Input:** a per-contrast `*_results.csv`  
**Output:** BrainNetViewer edge files (lead SNP per edge can be selected by smallest FDR/p)

```bash
bash run.sh python scripts/imggenetics/brainnetviewer_export.py   --results-csv outputs/imggenetics_post/GWAS_Results_by_Comparison/AsymAD_vs_Control_results.csv   --out-dir outputs/bnv   --keep-lowest-fdr-per-edge
```

---

## Stage 5 — blood differential expression (ADNI + ANMerge) (optional)

### ADNI blood expression preprocessing
**Input:** raw expression export + covariates table  
**Output:** cleaned subject-by-gene matrix + merged metadata

```bash
bash run.sh python scripts/transcriptomics/adni_expression_preprocess.py   --expression-file /path/to/ADNI_Gene_Expression_Profile.csv   --covariates-file /path/to/covariates.csv   --out-matrix outputs/transcriptomics/adni_expression_matrix_cleaned.csv   --out-metadata outputs/transcriptomics/adni_metadata_merged.csv
```

### ADNI DEG (covariate-adjusted)
```bash
bash run.sh python scripts/transcriptomics/adni_deg_covariates.py   --matrix outputs/transcriptomics/adni_expression_matrix_cleaned.csv   --metadata outputs/transcriptomics/adni_metadata_merged.csv   --out-dir outputs/transcriptomics/adni_deg
```

### ANMerge DEG (covariate-adjusted)
```bash
bash run.sh python scripts/transcriptomics/anmerge_deg_covariates.py   --input-csv /path/to/ANMerge_Unified_Ready_for_DEA.csv   --out-dir outputs/transcriptomics/anmerge_deg
```

---

## Stage 6 — ANMerge expression × MRI validation (optional)

### Prep merged expression + MRI tables for selected genes
```bash
bash run.sh python scripts/validation/anmerge_prep_expr_mri.py   --data-dir /path/to/anmerge_folder   --out-dir outputs/validation/results_b2
```

### Fit interaction models and apply within-gene FDR across regions
```bash
bash run.sh python scripts/validation/anmerge_mri_interaction_pergene_fdr.py   --data-dir outputs/validation/results_b2   --out-dir outputs/validation/results_b2_optionB_pergeneFDR   --qgene-threshold 0.05
```

---

## Notes on external tools

Some stages require external tools that must be installed separately:
- **FSL** (for rs-fMRI preprocessing)
- **PLINK** (for GWAS)
- **VEP** (variant consequence annotation; rsID lists are produced here)
- **Cytoscape + stringApp** (optional for PPI visualization)

---

## Acknowledgements

This code is designed for use with controlled-access datasets, including ADNI and ANMerge. If you use ADNI data, follow ADNI’s publication acknowledgement requirements. GTEx is referenced for optional cis-eQTL support in the broader project context.

---

## License

MIT (see `LICENSE`).
