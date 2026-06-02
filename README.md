# Subtype-Aware Imaging Genetics For rs-fMRI Connectivity

This repository contains the code used for a subtype-aware imaging-genetics workflow in Alzheimer's disease. The pipeline links GWAS-prioritized SNPs to resting-state fMRI connectivity edges, identifies candidate resilience-related connections, and supports downstream variant, transcriptomic, and MRI interpretation.

Protected ADNI and ANMerge datasets are not distributed with this repository. The code is written so every stage can be run from explicit input and output paths.

## Repository Layout

```text
.
├── src/imggenetics/                 # installable CLI for pairwise SNP-connectivity analysis
├── scripts/
│   ├── fmri/                        # fMRIPrep submission + MSDL connectivity extraction
│   ├── fmriphenotype/               # phenotype-only subtype model / candidate resilience edges
│   ├── gwas/                        # PLINK helpers and genotype-matrix construction
│   ├── imggenetics/                 # VEP and BrainNetViewer post-processing
│   ├── transcriptomics/             # ADNI + ANMerge blood-expression analyses
│   └── validation/                  # ANMerge expression x MRI interaction models
├── docs/                            # workflow notes and GWAS documentation
├── data/README.md                   # example local data layout
├── requirements.txt
├── pyproject.toml
└── run.sh                           # convenience wrapper for venv + execution
```

## Installation

```bash
git clone https://github.com/xiaoqinghuanglab/IG-fmri-snp.git
cd IG-fmri-snp

# recommended wrapper
bash run.sh python -c "import imggenetics; print('ok')"
```

Manual install:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

## Workflow Summary

1. `scripts/fmri/preprocess_fmri_fsl.sh`
   Runs participant-level `fMRIPrep 25.2.3` in a SLURM/Apptainer workflow.
2. `scripts/fmri/compute_connectivity_msdl.py`
   Extracts MSDL ROI time series from MNI-space fMRIPrep BOLD files, applies confound regression plus filtering, fits `GraphicalLassoCV`, and exports 741 signed partial-correlation edges per subject.
3. `scripts/fmriphenotype/get_resilience_edge.py`
   Fits phenotype-only subtype models to identify candidate resilience-related connectivity edges.
4. `scripts/gwas/*.py` and `scripts/gwas/run_gwas_plink.sh`
   Builds PLINK phenotype/covariate files, runs contrast-wise GWAS, extracts top SNPs, and constructs the union genotype matrix plus combined SNP metadata.
5. `imggenetics`
   Runs the pairwise SNP-connectivity interaction analysis for:
   - `AsymAD_vs_TypAD`
   - `TypAD_vs_Control`
   - `AsymAD_vs_Control`
6. `scripts/imggenetics/*.py`
   Exports VEP rsID lists and BrainNetViewer edge files from significant SNP-connectivity results.
7. `scripts/transcriptomics/` and `scripts/validation/`
   Optional ADNI / ANMerge follow-up analyses for gene-level interpretation.

## Standard Input Files

The main imaging-genetics stage expects these files:

| File | Required content |
|---|---|
| `msdl_all_subjects_connectivity_edges.csv` | `SubjectId` plus 741 `Connectivity_*` edge columns |
| `covariate_file.csv` | `SubjectId`, `Age`, `Sex`, `Scan_type`, `Manufacturer`, `Subtype` |
| `genotype_matrix.csv` | one row per subject plus SNP dosage columns |
| `snp_metadata.csv` | `SNP`, `Case`, `Control`, and optional `CHR`, `BP`, `P`, `OR` |
| `model2_candidate_resilience_edges.csv` | candidate edge file from the phenotype-only stage |

See [data/README.md](/Users/shihab/Codex/IG-fmri-snp/data/README.md:1) for an example local layout.

## Stage 1: fMRIPrep

The submission script in [scripts/fmri/preprocess_fmri_fsl.sh](/Users/shihab/Codex/IG-fmri-snp/scripts/fmri/preprocess_fmri_fsl.sh:1) is HPC-specific and should be edited for your dataset paths, FreeSurfer license path, and container location.

```bash
sbatch scripts/fmri/preprocess_fmri_fsl.sh
```

Expected outputs are BIDS derivatives under your fMRIPrep output directory, including MNI-space preprocessed resting-state BOLD files and matching confounds tables.

## Stage 2: MSDL Connectivity Extraction

This stage uses only the MSDL atlas and `GraphicalLassoCV`, matching the manuscript workflow.

```bash
bash run.sh python scripts/fmri/compute_connectivity_msdl.py \
  --input-dir /path/to/fmriprep_derivatives \
  --output-csv outputs/connectivity/msdl_all_subjects_connectivity_edges.csv
```

Key implementation details:

- input BOLD files are taken from `sub-*/func/*task-rest*_space-MNI152NLin2009cAsym*_desc-preproc_bold.nii.gz`
- nuisance regression uses common fMRIPrep motion and tissue confounds
- MSDL ROI time series are detrended, band-pass filtered, and z-scored
- connectivity is exported as signed partial-correlation edges from the `GraphicalLassoCV` precision matrix

Outputs:

- `msdl_all_subjects_connectivity_edges.csv`
- a processing log CSV next to the output matrix

## Stage 3: Phenotype-Only Resilience Edge Discovery

```bash
bash run.sh python scripts/fmriphenotype/get_resilience_edge.py \
  --connectivity outputs/connectivity/msdl_all_subjects_connectivity_edges.csv \
  --covariate /path/to/covariate_file.csv \
  --output_dir outputs/fmriphenotype
```

Important outputs:

- `model2_all_edge_results.csv.gz`
- `model2_candidate_resilience_edges.csv`
- `model2_strict_resilience_edges.csv`
- `model2_diagnostic_report.txt`

## Stage 4: GWAS Feature Prioritization

1. Build per-contrast phenotype and covariate files:

```bash
bash run.sh python scripts/gwas/prepare_pheno_covar.py \
  --clinical /path/to/clinical_covariates.csv \
  --out-dir outputs/gwas/pheno_covar
```

2. Run PLINK logistic regression:

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

Repeat for:

- `AsymAD_vs_TypAD`
- `TypAD_vs_Control`

4. Build the union genotype matrix and combined SNP metadata:

```bash
bash run.sh python scripts/gwas/make_genotype_matrix.py \
  --bfile /path/to/QC_GENOTYPES_PREFIX \
  --top-list outputs/gwas/top1000_AsymAD_vs_Control.tsv AsymAD Control \
  --top-list outputs/gwas/top1000_AsymAD_vs_TypAD.tsv AsymAD TypAD \
  --top-list outputs/gwas/top1000_TypAD_vs_Control.tsv TypAD Control \
  --out-dir outputs/gwas/final
```

Outputs:

- `genotype_matrix.csv`
- `snp_metadata.csv`

## Stage 5: Pairwise SNP-Connectivity Analysis

The `imggenetics` CLI now runs the pairwise interaction workflow used in the paper.

```bash
bash run.sh imggenetics \
  --connectivity-file outputs/connectivity/msdl_all_subjects_connectivity_edges.csv \
  --covariates-file /path/to/covariate_file.csv \
  --genetic-file outputs/gwas/final/genotype_matrix.csv \
  --snp-metadata-file outputs/gwas/final/snp_metadata.csv \
  --candidate-edges-file outputs/fmriphenotype/model2_candidate_resilience_edges.csv \
  --output-dir outputs/imggenetics
```

Per comparison, the CLI writes:

- `<comparison>_full_results.csv.gz`
- `<comparison>_significant_interaction_fdr05.csv`
- `<comparison>_significant_interaction_fdr05_minimal.csv`
- `<comparison>_lead_snp_per_edge_fdr05.csv`
- `<comparison>_summary_counts.csv`
- `<comparison>_skipped_snps.csv`
- `comparison_report.txt`

Comparison folders:

- `AsymAD_vs_TypAD/`
- `TypAD_vs_Control/`
- `AsymAD_vs_Control/`

## Stage 6: Post-Processing

### VEP Input Lists

```bash
bash run.sh python scripts/imggenetics/export_vep_inputs.py \
  --results-dir outputs/imggenetics \
  --out-dir outputs/imggenetics_post
```

Outputs:

- `VEP_Input_Files/<comparison>_VEP.txt`
- `GWAS_Results_by_Comparison/<comparison>_results.csv`

### BrainNetViewer Edge Files

Use either the full significant table or the lead-SNP-per-edge table for one comparison.

```bash
bash run.sh python scripts/imggenetics/brainnetviewer_export.py \
  --results-csv outputs/imggenetics/AsymAD_vs_TypAD/AsymAD_vs_TypAD_lead_snp_per_edge_fdr05.csv \
  --out-dir outputs/brainnet \
  --label AsymAD_vs_TypAD
```

Outputs:

- `msdl_nodes.node`
- `edges_<label>.edge`

## Optional Follow-Up Analyses

- [scripts/transcriptomics/adni_expression_preprocess.py](/Users/shihab/Codex/IG-fmri-snp/scripts/transcriptomics/adni_expression_preprocess.py:1)
- [scripts/transcriptomics/adni_deg_covariates.py](/Users/shihab/Codex/IG-fmri-snp/scripts/transcriptomics/adni_deg_covariates.py:1)
- [scripts/transcriptomics/anmerge_deg_covariates.py](/Users/shihab/Codex/IG-fmri-snp/scripts/transcriptomics/anmerge_deg_covariates.py:1)
- [scripts/validation/anmerge_prep_expr_mri.py](/Users/shihab/Codex/IG-fmri-snp/scripts/validation/anmerge_prep_expr_mri.py:1)
- [scripts/validation/anmerge_mri_interaction_pergene_fdr.py](/Users/shihab/Codex/IG-fmri-snp/scripts/validation/anmerge_mri_interaction_pergene_fdr.py:1)

## Notes

- This repository is intentionally code-only; notebooks and legacy scripts were removed to keep the paper submission workflow clean.
- ADNI / ANMerge data remain external and must be supplied locally.
- The paired imaging-genetics sample is modest, especially for Typical AD, so outputs should be interpreted as hypothesis-generating unless independently replicated.
