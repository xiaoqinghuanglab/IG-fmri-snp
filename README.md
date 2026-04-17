# Subtype-aware imaging genetics for rs-fMRI connectivity

This repository implements a subtype-aware imaging-genetics pipeline for Alzheimer’s disease, linking common genetic variation to resting-state fMRI functional connectivity. It harmonizes ADNI genotypes, prioritizes SNPs via GWAS, derives MSDL-based connectivity features, runs phenotype-only subtype analysis to identify resilience-related edges, and fits SNP×subtype OLS models with FDR control. Downstream modules support variant annotation, enrichment, PPI, and cross-cohort validation.

This repository collects a complete **imaging-genetics workflow** that connects:

- **Genetics:** GWAS-based SNP feature prioritization per subtype contrast (PLINK)
- **Imaging:** rs-fMRI preprocessing with **fMRIPrep** and connectivity feature extraction (MSDL 39 ROI → 741 edges)
- **Phenotype analysis:** subtype-only edge modeling to identify **candidate and strict resilience edges**
- **Association testing:** SNP × subtype interaction models for SNP–edge pairs (Python / statsmodels)
- **Downstream interpretation:** VEP rsID lists, BrainNetViewer edge exports, enrichment/PPI helpers
- **External validation (optional):** blood transcriptomics (ADNI + ANMerge) and ANMerge expression × MRI models

All runnable steps are exposed as **command-line scripts** with **explicit input/output path flags** so users can run the same code without editing core source files.

> **Protected data:** ADNI/ANMerge datasets are not distributed with this repo. See `data/README.md` for an example local layout and required filenames/columns.

---

## Repository layout

```text
.
├── src/imggenetics/                 # installable package: SNP×Subtype OLS for connectivity phenotypes
├── scripts/
│   ├── fmri/                        # fMRIPrep preprocessing + MSDL connectivity derivation
│   ├── fmriphenotype/               # phenotype-only subtype analysis + resilience edge discovery
│   │   └── get_resilience_edge.py
│   ├── gwas/                        # PLINK helpers + top-SNP extraction + genotype matrix builder
│   ├── imggenetics/                 # post-processing utilities (VEP lists, BrainNetViewer exports)
│   ├── transcriptomics/             # ADNI + ANMerge blood DEG analysis
│   └── validation/                  # ANMerge expression × MRI interaction models
├── notebooks/                       # original notebook version (reference)
├── docs/                            # additional workflow notes + GWAS documentation PDF
├── data/README.md                   # data layout (local only)
├── requirements.txt                 # Python dependencies
└── run.sh                           # convenience runner (creates .venv, installs deps, runs commands)

# Subtype-aware imaging genetics for rs-fMRI connectivity

This repository implements a subtype-aware imaging-genetics pipeline for Alzheimer’s disease, linking common genetic variation to resting-state fMRI functional connectivity. It harmonizes ADNI genotypes, prioritizes SNPs via GWAS, derives MSDL-based connectivity features, runs phenotype-only subtype analysis to identify resilience-related edges, and fits SNP×subtype OLS models with FDR control. Downstream modules support variant annotation, enrichment, PPI, and cross-cohort validation.

This repository collects a complete **imaging-genetics workflow** that connects:

- **Genetics:** GWAS-based SNP feature prioritization per subtype contrast (PLINK)
- **Imaging:** rs-fMRI preprocessing with **fMRIPrep** and connectivity feature extraction (MSDL 39 ROI → 741 edges)
- **Phenotype analysis:** subtype-only edge modeling to identify **candidate and strict resilience edges**
- **Association testing:** SNP × subtype interaction models for SNP–edge pairs (Python / statsmodels)
- **Downstream interpretation:** VEP rsID lists, BrainNetViewer edge exports, enrichment/PPI helpers
- **External validation (optional):** blood transcriptomics (ADNI + ANMerge) and ANMerge expression × MRI models

All runnable steps are exposed as **command-line scripts** with **explicit input/output path flags** so users can run the same code without editing core source files.

> **Protected data:** ADNI/ANMerge datasets are not distributed with this repo. See `data/README.md` for an example local layout and required filenames/columns.

---

## Repository layout

```text
.
├── src/imggenetics/                 # installable package: SNP×Subtype OLS for connectivity phenotypes
├── scripts/
│   ├── fmri/                        # fMRIPrep preprocessing + MSDL connectivity derivation
│   ├── fmriphenotype/               # phenotype-only subtype analysis + resilience edge discovery
│   │   └── get_resilience_edge.py
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
- an **edge feature matrix** (subjects × connectivity edges),
- **subject covariates** (age/sex/scanner/scan_type + subtype label),
- and a **SNP metadata table** that groups SNPs by their originating GWAS contrast,

this stage fits an OLS model per SNP–edge pair with a genotype-by-subtype interaction and writes:

- full model results
- BH-FDR corrected p-values
- a filtered significant table

### Required inputs

You can place these 4 CSVs anywhere and pass explicit paths:

| File | Required columns (minimum) |
|---|---|
| `genotype_matrix.csv` | `Subject_ID` + SNP columns (for example, `rs12345_A`) |
| `Connectivity_matrix_all_subjects_region_pairs.csv` | `Subject_ID` + connectivity columns |
| `subject_covariates.csv` | `Subject_ID`, `Subtype`, `age`, `sex`, `scanner`, `scan_type` |
| `snp_metadata.csv` | `SNP` (rsID), `Case`, `Control` |

### Run

```bash
bash run.sh imggenetics \
  --output-dir outputs/imggenetics \
  --genetic-file /path/to/genotype_matrix.csv \
  --connectivity-file /path/to/Connectivity_matrix_all_subjects_region_pairs.csv \
  --covariates-file /path/to/subject_covariates.csv \
  --snp-metadata-file /path/to/snp_metadata.csv \
  --subtypes Control TypAD AsymAD \
  --reference-subtype Control
```

### Outputs

Written under `--output-dir`:

- `all_regression_result.csv`
- `all_regression_result_FDR.csv`
- `significant_result.csv`

---

## Stage 1 — rs-fMRI preprocessing with fMRIPrep

### What this step does

Resting-state fMRI preprocessing is performed with **fMRIPrep** in a SLURM array-job workflow. The preprocessing script is designed for BIDS-formatted input data and runs each participant independently through an Apptainer container.

The workflow includes:

- participant-wise preprocessing with **fMRIPrep**
- FreeSurfer template setup (`fsaverage`, `fsaverage5`) before parallel jobs proceed
- BIDS-aware derivatives output
- support for both:
  - **full functional preprocessing** when BOLD data are available
  - **anatomical-only processing** when a subject has no BOLD file

### Key implementation details

The current preprocessing workflow:

- uses **fMRIPrep 25.2.3**
- expects a **BIDS dataset**
- writes outputs to a derivatives directory
- uses **MNI152NLin2009cAsym:res-2**, `anat`, `func`, and `fsaverage5` output spaces
- supports `--cifti-output`
- uses fieldmap-less susceptibility distortion correction with `--use-syn-sdc warn`
- is configured for HPC execution with **SLURM + Apptainer**

### Requirements

Before running this stage, make sure you have:

- a valid **BIDS dataset**
- a **FreeSurfer license file**
- an **fMRIPrep container image**
- SLURM access if running on HPC

### Typical inputs

- `ADNI_Outputs/datasets/<DATASET_NAME>/`
- `participants.tsv` inside that BIDS dataset
- subject folders such as `sub-XXXX/...`

### Typical outputs

Written under:

- `ADNI_Outputs/derivatives/<DATASET_NAME>/`

These derivatives are then used for downstream connectivity extraction.

### Run

Edit the dataset- and environment-specific paths in the SLURM script under `scripts/fmri/`, then submit it with:

```bash
sbatch scripts/fmri/<your_fmriprep_slurm_script>.sh
```

### Notes

- If a subject has no BOLD file, the script automatically switches to **anat-only** mode.
- If FreeSurfer output already exists for a subject, that subject is skipped.
- The script is designed for large cohort processing using a job array.

---

## Stage 2 — connectivity features (MSDL + GraphicalLassoCV)

### What this step does

For each subject’s preprocessed fMRI file:

- extracts ROI time series using the **MSDL atlas (39 ROIs)**
- optionally aligns time-series length (truncate/pad)
- fits **GraphicalLassoCV**
- exports edge-wise connectivity features

This produces **741** edge features per subject (39×38/2).

### Inputs

- preprocessed fMRI derivatives from Stage 1
- filenames containing the subject identifier, or a user-provided regex for ID extraction

### Outputs

- a subject-by-edge CSV with one row per subject and one column per connectivity edge

### Run

```bash
bash run.sh python scripts/fmri/compute_connectivity_msdl.py \
  --input-dir /path/to/preprocessed_fmri \
  --pattern "*_bold.nii.gz" \
  --output-csv outputs/connectivity/Connectivity_matrix_all_subjects_region_pairs.csv \
  --tr 3 \
  --fixed-length median \
  --abs-edges
```

---

## Stage 3 — phenotype-only subtype analysis and resilience edge discovery

### What this step does

This stage runs a **phenotype-only subtype model** on connectivity edges before the SNP-level imaging-genetics association step. Its purpose is to identify edges whose group-level patterns are consistent with:

- **resilience**
- **compensation**
- **non-resilience / pathology-like patterns**
- **mixed or unclear patterns**

This step is especially useful for prioritizing connectivity edges for resilience-focused downstream interpretation.

### Model

For each edge, the script fits a covariate-adjusted linear model of the form:

```text
edge ~ C(Subtype) + Age + Sex + C(Scan_type) + C(Manufacturer)
```

with **Control** as the reference group.

### What the script evaluates

For each edge, the script computes:

- omnibus subtype test
- `AsymAD vs Control`
- `TypicalAD vs Control`
- `AsymAD vs TypicalAD`
- FDR-corrected p-values
- adjusted subtype means
- pattern labels based on group ordering and significance

### Pattern categories

The script produces both **strict** and **candidate** pattern labels, including:

- **Resilience**
  - `AsymAD ≈ Control > TypicalAD`
  - `AsymAD ≈ Control < TypicalAD`
- **Compensation**
  - AsymAD differs from both other groups
- **Non-resilience**
  - `AsymAD ≈ TypicalAD`, both different from Control
- **Mixed / unclear**

### Required inputs

| File | Required columns (minimum) |
|---|---|
| `msdl_all_subjects_connectivity_edges.csv` | subject ID column + edge columns |
| `covariate_file.csv` | `SubjectId`, `Age`, `Sex`, `Scan_type`, `Manufacturer`, `Subtype` |

### Run

```bash
bash run.sh python scripts/fmriphenotype/get_resilience_edge.py \
  --connectivity /path/to/msdl_all_subjects_connectivity_edges.csv \
  --covariate /path/to/covariate_file.csv \
  --output_dir outputs/fmriphenotype
```

### Outputs

Written under `--output_dir`:

- `model2_all_edge_results.csv.gz`
- `model2_strict_resilience_edges.csv`
- `model2_candidate_resilience_edges.csv`
- `model2_strict_compensation_edges.csv`
- `model2_candidate_compensation_edges.csv`
- `model2_strict_non_resilience_edges.csv`
- `model2_candidate_non_resilience_edges.csv`
- `model2_strict_mixed_unclear_edges.csv`
- `model2_candidate_mixed_unclear_edges.csv`
- `model2_skipped_edges.csv`
- `model2_top50_by_omnibus_fdr.csv`
- `model2_top50_by_raw_typ_vs_control.csv`
- `model2_diagnostic_report.txt`

### Why this stage is useful

This phenotype-based screening layer helps identify **candidate resilience edges** before SNP-level interpretation. It provides a biologically motivated way to focus on edges where AsymAD remains closer to Control while TypicalAD diverges.

---

## Stage 4 — GWAS for SNP feature prioritization (PLINK)

### What this stage is for

In this project, GWAS is used as a **feature selection step**: per contrast, SNPs are ranked by association p-value and the top set is carried forward into the SNP×connectivity modeling.

A full, command-oriented GWAS walkthrough is included in:

- `docs/GWAS_Analysis_Documentation.pdf`
- `docs/gwas_workflow.md`

### Typical inputs

- QC-passed PLINK files: `<PREFIX>.bed/.bim/.fam`
- a clinical/covariate table with subject IDs, subtype labels, age, sex, APOE4, and ancestry PCs

### Typical outputs

- contrast-specific PLINK association results (`*.assoc.logistic`)
- `top1000_<contrast>.tsv`
- `snp_metadata.csv`
- `genotype_matrix.csv`

### Minimal run

```bash
# 1) build per-contrast phenotype + covariate files
bash run.sh python scripts/gwas/prepare_pheno_covar.py \
  --clinical-csv /path/to/clinical_covariates.csv \
  --out-dir outputs/gwas/pheno_covar

# 2) run PLINK logistic regression per contrast
bash run.sh bash scripts/gwas/run_gwas_plink.sh \
  --bfile /path/to/QC_GENOTYPES_PREFIX \
  --pheno-dir outputs/gwas/pheno_covar \
  --out-dir outputs/gwas/assoc

# 3) extract top SNPs per contrast
bash run.sh python scripts/gwas/extract_top_snps.py \
  --assoc-dir outputs/gwas/assoc \
  --out-dir outputs/gwas \
  --top-n 1000

# 4) build genotype_matrix.csv + snp_metadata.csv
bash run.sh python scripts/gwas/make_genotype_matrix.py \
  --bfile /path/to/QC_GENOTYPES_PREFIX \
  --top-list outputs/gwas/top1000_AsymAD_vs_Control.tsv AsymAD Control \
  --top-list outputs/gwas/top1000_AsymAD_vs_TypAD.tsv AsymAD TypAD \
  --top-list outputs/gwas/top1000_TypAD_vs_Control.tsv TypAD Control \
  --out-dir outputs/gwas/final
```

---

## Stage 5 — post-processing (VEP lists + BrainNetViewer edge files)

These are optional utilities that turn `imggenetics` outputs into downstream-friendly files.

### VEP rsID lists + per-contrast result CSVs

**Input:** `significant_result.csv`  
**Output:** `VEP_Input_Files/*.txt` + `GWAS_Results_by_Comparison/*_results.csv`

```bash
bash run.sh python scripts/imggenetics/export_vep_inputs.py \
  --significant-csv outputs/imggenetics/significant_result.csv \
  --out-dir outputs/imggenetics_post
```

### BrainNetViewer edge exports

**Input:** a per-contrast `*_results.csv`  
**Output:** BrainNetViewer edge files

```bash
bash run.sh python scripts/imggenetics/brainnetviewer_export.py \
  --results-csv outputs/imggenetics_post/GWAS_Results_by_Comparison/AsymAD_vs_Control_results.csv \
  --out-dir outputs/bnv \
  --keep-lowest-fdr-per-edge
```

---

## Stage 6 — blood differential expression (ADNI + ANMerge) (optional)

### ADNI blood expression preprocessing

```bash
bash run.sh python scripts/transcriptomics/adni_expression_preprocess.py \
  --expression-file /path/to/ADNI_Gene_Expression_Profile.csv \
  --covariates-file /path/to/covariates.csv \
  --out-matrix outputs/transcriptomics/adni_expression_matrix_cleaned.csv \
  --out-metadata outputs/transcriptomics/adni_metadata_merged.csv
```

### ADNI DEG (covariate-adjusted)

```bash
bash run.sh python scripts/transcriptomics/adni_deg_covariates.py \
  --matrix outputs/transcriptomics/adni_expression_matrix_cleaned.csv \
  --metadata outputs/transcriptomics/adni_metadata_merged.csv \
  --out-dir outputs/transcriptomics/adni_deg
```

### ANMerge DEG (covariate-adjusted)

```bash
bash run.sh python scripts/transcriptomics/anmerge_deg_covariates.py \
  --input-csv /path/to/ANMerge_Unified_Ready_for_DEA.csv \
  --out-dir outputs/transcriptomics/anmerge_deg
```

---

## Stage 7 — ANMerge expression × MRI validation (optional)

### Prep merged expression + MRI tables for selected genes

```bash
bash run.sh python scripts/validation/anmerge_prep_expr_mri.py \
  --data-dir /path/to/anmerge_folder \
  --out-dir outputs/validation/results_b2
```

### Fit interaction models and apply within-gene FDR across regions

```bash
bash run.sh python scripts/validation/anmerge_mri_interaction_pergene_fdr.py \
  --data-dir outputs/validation/results_b2 \
  --out-dir outputs/validation/results_b2_optionB_pergeneFDR \
  --qgene-threshold 0.05
```

---

## Notes on external tools

Some stages require external tools that must be installed separately:

- **fMRIPrep** for rs-fMRI preprocessing
- **Apptainer / Singularity** for containerized HPC execution
- **FreeSurfer license** for fMRIPrep structural processing
- **PLINK** for GWAS
- **VEP** for variant consequence annotation
- **Cytoscape + stringApp** (optional) for PPI visualization

---

## Acknowledgements

This code is designed for use with controlled-access datasets, including ADNI and ANMerge. If you use ADNI data, follow ADNI’s publication acknowledgement requirements. GTEx is referenced for optional cis-eQTL support in the broader project context.

---

## License

MIT (see `LICENSE`).
