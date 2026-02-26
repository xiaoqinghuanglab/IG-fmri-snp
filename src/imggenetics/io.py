"""I/O helpers for the imggenetics pipeline.

This module is intentionally small and boring:
  - load the four required CSV inputs (either from a data directory or
    from explicit paths)
  - standardize the subject identifier column name to `Subject_ID`

Keeping I/O logic in one place makes the rest of the pipeline easier to test
and reduces "path editing" for users.
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import pandas as pd


def _read_csv(path: Path) -> pd.DataFrame:
    """Thin wrapper around pandas.read_csv.

    Having a single function makes it easy to:
      - add shared read_csv options later (dtypes, NA values, etc.)
      - swap to parquet/feather in the future without touching the pipeline
    """
    return pd.read_csv(path)


def load_inputs_from_paths(
    genetic_path: Path,
    connectivity_path: Path,
    covariates_path: Path,
    snp_metadata_path: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load the 4 required CSV inputs from explicit file paths."""
    g = _read_csv(genetic_path)
    c = _read_csv(connectivity_path)
    cov = _read_csv(covariates_path)
    meta = _read_csv(snp_metadata_path)
    return g, c, cov, meta


def load_inputs(
    data_dir: Path,
    genetic_name: str = "genotype_matrix.csv",
    connectivity_name: str = "Connectivity_matrix_all_subjects_region_pairs.csv",
    covariates_name: str = "subject_covariates.csv",
    snp_metadata_name: str = "snp_metadata.csv",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load the 4 required CSV inputs from a folder using default filenames."""
    g = _read_csv(data_dir / genetic_name)
    c = _read_csv(data_dir / connectivity_name)
    cov = _read_csv(data_dir / covariates_name)
    meta = _read_csv(data_dir / snp_metadata_name)
    return g, c, cov, meta


def harmonize_subject_id(
    genetic: pd.DataFrame,
    connectivity: pd.DataFrame,
    covariates: pd.DataFrame,
) -> None:
    """Harmonize subject ID column names to the shared name: `Subject_ID`.

    Different exports often use slightly different headers. We normalize a few
    common variants in-place:
      - connectivity: 'Subject ID' -> 'Subject_ID'
      - covariates:   'subject_id' -> 'Subject_ID'
      - genetic:      'SubjectId'  -> 'Subject_ID'

    The remainder of the pipeline assumes the join key is exactly `Subject_ID`.
    """

    # Rename common variants (in-place).
    if "Subject ID" in connectivity.columns:
        connectivity.rename(columns={"Subject ID": "Subject_ID"}, inplace=True)
    if "subject_id" in covariates.columns:
        covariates.rename(columns={"subject_id": "Subject_ID"}, inplace=True)
    if "SubjectId" in genetic.columns:
        genetic.rename(columns={"SubjectId": "Subject_ID"}, inplace=True)

    # Validate all inputs have the required join column.
    for name, df in [
        ("genetic", genetic),
        ("connectivity", connectivity),
        ("covariates", covariates),
    ]:
        if "Subject_ID" not in df.columns:
            raise KeyError(
                f"Missing 'Subject_ID' in {name} dataframe after harmonization. "
                f"Found columns: {list(df.columns)[:20]}"
            )
