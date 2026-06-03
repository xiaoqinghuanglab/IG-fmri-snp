# fMRI connectivity workflow

Use the scripts in this folder when you want to regenerate the connectivity input from raw fMRIPrep outputs instead of using the bundled `data/msdl_all_subjects_connectivity_edges.csv`.

## Files

- `preprocess_fmri_fsl.sh`
  HPC-oriented SLURM / Apptainer launcher for `fMRIPrep 25.2.3`.
- `compute_connectivity_msdl.py`
  MSDL-only connectivity extraction with `GraphicalLassoCV`, matching the paper workflow.

## Regeneration path

1. Run `fMRIPrep` on resting-state BOLD data.
2. Point `compute_connectivity_msdl.py` at the derivatives directory.

Example:

```bash
bash run.sh python scripts/fmri/compute_connectivity_msdl.py \
  --input-dir /path/to/fmriprep_derivatives \
  --output-csv outputs/connectivity/msdl_all_subjects_connectivity_edges.csv
```

The exported CSV should contain one row per subject and 741 signed MSDL partial-correlation edges.

## Notes

- The repo's cleaned downstream workflow uses only the MSDL output, not alternate atlases.
- `preprocess_fmri_fsl.sh` keeps its historical filename, but it launches `fMRIPrep`, not an FSL-only preprocessing pipeline.
- You do not need to run this folder to reproduce the paper results bundled with the repo; the final MSDL connectivity CSV is already included in `data/`.
