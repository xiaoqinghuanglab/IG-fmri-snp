from __future__ import annotations

import pandas as pd
from pathlib import Path
from typing import Tuple


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def load_inputs(
    data_dir: Path,
    genetic_name: str = "genotype_matrix.csv",
    connectivity_name: str = "Connectivity_matrix_all_subjects_region_pairs.csv",
    covariates_name: str = "subject_covariates.csv",
    snp_metadata_name: str = "snp_metadata.csv",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load the 4 required CSV inputs.
    """
    g = _read_csv(data_dir / genetic_name)
    c = _read_csv(data_dir / connectivity_name)
    cov = _read_csv(data_dir / covariates_name)
    meta = _read_csv(data_dir / snp_metadata_name)
    return g, c, cov, meta


def harmonize_subject_id(genetic: pd.DataFrame, connectivity: pd.DataFrame, covariates: pd.DataFrame):
    """
    Harmonize subject ID column names to 'Subject_ID' to match the notebook.

    Supported variants:
      - connectivity: 'Subject ID' -> 'Subject_ID'
      - covariates: 'subject_id' -> 'Subject_ID'
      - genetic: 'SubjectId' -> 'Subject_ID'
    """
    if "Subject ID" in connectivity.columns:
        connectivity.rename(columns={"Subject ID": "Subject_ID"}, inplace=True)
    if "subject_id" in covariates.columns:
        covariates.rename(columns={"subject_id": "Subject_ID"}, inplace=True)
    if "SubjectId" in genetic.columns:
        genetic.rename(columns={"SubjectId": "Subject_ID"}, inplace=True)

    # If still missing, raise a clear error
    for name, df in [("genetic", genetic), ("connectivity", connectivity), ("covariates", covariates)]:
        if "Subject_ID" not in df.columns:
            raise KeyError(f"Missing 'Subject_ID' in {name} dataframe after harmonization. "
                           f"Found columns: {list(df.columns)[:20]}")
