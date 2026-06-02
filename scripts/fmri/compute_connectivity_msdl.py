#!/usr/bin/env python3
"""Extract MSDL GraphicalLassoCV connectivity edges from fMRIPrep outputs.

This script follows the manuscript-aligned connectivity workflow:
  - uses fMRIPrep MNI-space preprocessed resting-state BOLD files
  - applies nuisance regression using common fMRIPrep confounds
  - uses the 39-region MSDL atlas
  - extracts ROI time series with detrending, filtering, and z-scoring
  - estimates GraphicalLassoCV connectivity
  - converts the precision matrix to signed partial correlations
  - exports one subject-by-edge CSV with 741 upper-triangle edges
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from nilearn.datasets import fetch_atlas_msdl
from nilearn.maskers import NiftiMapsMasker
from sklearn.covariance import GraphicalLassoCV
from sklearn.exceptions import ConvergenceWarning


MSDL_EXPECTED_ROIS = 39


def normalize_subject_id(value: str) -> str:
    text = str(value).strip()
    if text.startswith("sub-"):
        text = text[4:]

    if re.match(r"^\d{3}_S_\d{4}$", text):
        return text

    match = re.match(r"^(\d{3})S(\d{4})$", text)
    if match:
        return f"{match.group(1)}_S_{match.group(2)}"

    return text


def parse_args():
    parser = argparse.ArgumentParser(description="Compute MSDL partial-correlation connectivity from fMRIPrep outputs.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="fMRIPrep derivatives directory containing sub-*/func folders.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        required=True,
        help="Output CSV path for the subject-by-edge matrix.",
    )
    parser.add_argument(
        "--log-csv",
        type=Path,
        default=None,
        help="Optional processing log CSV. Defaults to <output>.log.csv.",
    )
    parser.add_argument(
        "--atlas-cache-dir",
        type=Path,
        default=None,
        help="Optional cache directory for nilearn atlas downloads.",
    )
    parser.add_argument("--low-pass", type=float, default=0.10, help="Low-pass cutoff in Hz.")
    parser.add_argument("--high-pass", type=float, default=0.01, help="High-pass cutoff in Hz.")
    parser.add_argument(
        "--include-global-signal",
        action="store_true",
        help="Include global signal in the confound regression set.",
    )
    parser.add_argument(
        "--abs-edges",
        action="store_true",
        help="Take absolute values of partial-correlation edges before export.",
    )
    return parser.parse_args()


def read_json_value(json_file: Path, key: str, default=None):
    if not json_file.exists():
        return default
    try:
        with json_file.open("r") as handle:
            payload = json.load(handle)
        return payload.get(key, default)
    except Exception:
        return default


def get_subject_dirs(base_dir: Path) -> list[Path]:
    return sorted(path for path in base_dir.glob("sub-*") if path.is_dir())


def get_run_prefix_from_bold_name(bold_name: str) -> str:
    return bold_name.split("_space-")[0]


def find_bold_files(func_dir: Path) -> list[Path]:
    pattern = "*task-rest*_space-MNI152NLin2009cAsym*_desc-preproc_bold.nii.gz"
    return sorted(func_dir.glob(pattern))


def load_confounds(confounds_tsv: Path, include_global_signal: bool = False) -> pd.DataFrame | None:
    if not confounds_tsv.exists():
        return None

    df = pd.read_csv(confounds_tsv, sep="\t")
    base_cols = [
        "trans_x",
        "trans_y",
        "trans_z",
        "rot_x",
        "rot_y",
        "rot_z",
        "trans_x_derivative1",
        "trans_y_derivative1",
        "trans_z_derivative1",
        "rot_x_derivative1",
        "rot_y_derivative1",
        "rot_z_derivative1",
        "white_matter",
        "csf",
        "framewise_displacement",
        "std_dvars",
    ]
    if include_global_signal:
        base_cols.append("global_signal")

    selected = [column for column in base_cols if column in df.columns]
    selected.extend(column for column in df.columns if column.startswith("cosine"))
    selected.extend(column for column in df.columns if column.startswith("non_steady_state_outlier"))

    if not selected:
        return None

    confounds = df[selected].copy()
    confounds = confounds.fillna(0.0)
    return confounds


def precision_to_partial_corr(precision: np.ndarray) -> np.ndarray:
    diagonal = np.sqrt(np.diag(precision))
    outer = np.outer(diagonal, diagonal)
    with np.errstate(divide="ignore", invalid="ignore"):
        partial = -precision / outer
    partial[np.isnan(partial)] = 0.0
    partial[np.isinf(partial)] = 0.0
    np.fill_diagonal(partial, 1.0)
    return partial


def build_edge_features(partial_corr: np.ndarray, labels: list[str], abs_edges: bool) -> dict[str, float]:
    features: dict[str, float] = {}
    n_rois = partial_corr.shape[0]
    for i in range(n_rois):
        for j in range(i + 1, n_rois):
            value = float(partial_corr[i, j])
            if abs_edges:
                value = float(abs(value))
            features[f"Connectivity_{labels[i]}_{labels[j]}"] = value
    return features


def extract_run_features(
    bold_file: Path,
    confounds_file: Path,
    atlas,
    low_pass: float,
    high_pass: float,
    include_global_signal: bool,
    abs_edges: bool,
) -> dict[str, float]:
    run_prefix = get_run_prefix_from_bold_name(bold_file.name)
    bold_json = bold_file.with_suffix("").with_suffix(".json")
    tr = read_json_value(bold_json, "RepetitionTime", None)
    confounds = load_confounds(confounds_file, include_global_signal=include_global_signal)

    masker = NiftiMapsMasker(
        maps_img=atlas.maps,
        standardize="zscore_sample",
        detrend=True,
        low_pass=low_pass if tr is not None else None,
        high_pass=high_pass if tr is not None else None,
        t_r=tr,
        verbose=0,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)
        timeseries = masker.fit_transform(str(bold_file), confounds=confounds)

    labels = [str(label) for label in atlas.labels]
    if timeseries.shape[1] != MSDL_EXPECTED_ROIS or len(labels) != MSDL_EXPECTED_ROIS:
        raise ValueError(
            f"MSDL ROI count mismatch for {run_prefix}: got {timeseries.shape[1]}, expected {MSDL_EXPECTED_ROIS}."
        )

    model = GraphicalLassoCV(cv=5, max_iter=2000, n_jobs=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)
        warnings.simplefilter("ignore", category=UserWarning)
        model.fit(timeseries)

    partial_corr = precision_to_partial_corr(model.precision_)
    return build_edge_features(partial_corr, labels, abs_edges=abs_edges)


def main():
    args = parse_args()

    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    input_dir = args.input_dir.expanduser().resolve()
    output_csv = args.output_csv.expanduser().resolve()
    log_csv = args.log_csv.expanduser().resolve() if args.log_csv else output_csv.with_suffix(".log.csv")

    if not input_dir.exists():
        raise SystemExit(f"Input directory not found: {input_dir}")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    log_csv.parent.mkdir(parents=True, exist_ok=True)

    atlas = fetch_atlas_msdl(
        data_dir=str(args.atlas_cache_dir.expanduser().resolve()) if args.atlas_cache_dir else None,
        verbose=0,
    )
    subjects = get_subject_dirs(input_dir)
    if not subjects:
        raise SystemExit(f"No subject folders found under: {input_dir}")

    rows: list[dict[str, float | str]] = []
    log_rows: list[dict[str, object]] = []

    print(f"Found {len(subjects)} subject folders in {input_dir}")

    for subject_dir in subjects:
        func_dir = subject_dir / "func"
        subject_id = normalize_subject_id(subject_dir.name)

        if not func_dir.exists():
            log_rows.append(
                {
                    "subject_id": subject_id,
                    "status": "skip_no_func_folder",
                    "n_runs_found": 0,
                    "n_runs_used": 0,
                    "message": "No func folder",
                }
            )
            continue

        bold_files = find_bold_files(func_dir)
        if not bold_files:
            log_rows.append(
                {
                    "subject_id": subject_id,
                    "status": "skip_no_mni_bold",
                    "n_runs_found": 0,
                    "n_runs_used": 0,
                    "message": "No MNI-space preprocessed BOLD files found",
                }
            )
            continue

        run_features: list[dict[str, float]] = []
        run_errors: list[str] = []

        for bold_file in bold_files:
            run_prefix = get_run_prefix_from_bold_name(bold_file.name)
            confounds_file = func_dir / f"{run_prefix}_desc-confounds_timeseries.tsv"

            try:
                if not confounds_file.exists():
                    raise FileNotFoundError(f"Missing confounds file: {confounds_file.name}")

                features = extract_run_features(
                    bold_file=bold_file,
                    confounds_file=confounds_file,
                    atlas=atlas,
                    low_pass=args.low_pass,
                    high_pass=args.high_pass,
                    include_global_signal=args.include_global_signal,
                    abs_edges=args.abs_edges,
                )
                run_features.append(features)
            except Exception as exc:
                run_errors.append(f"{run_prefix}: {exc}")

        if not run_features:
            log_rows.append(
                {
                    "subject_id": subject_id,
                    "status": "error",
                    "n_runs_found": len(bold_files),
                    "n_runs_used": 0,
                    "message": " | ".join(run_errors),
                }
            )
            continue

        run_df = pd.DataFrame(run_features)
        subject_row = {"SubjectId": subject_id}
        for column in run_df.columns:
            subject_row[column] = float(run_df[column].mean())
        rows.append(subject_row)

        log_rows.append(
            {
                "subject_id": subject_id,
                "status": "success" if not run_errors else "partial_success",
                "n_runs_found": len(bold_files),
                "n_runs_used": len(run_features),
                "message": " | ".join(run_errors),
            }
        )

    if not rows:
        raise SystemExit("No subjects produced usable connectivity features.")

    df = pd.DataFrame(rows)
    edge_cols = sorted(column for column in df.columns if column != "SubjectId")
    df = df[["SubjectId"] + edge_cols].sort_values("SubjectId").reset_index(drop=True)
    df.to_csv(output_csv, index=False)
    pd.DataFrame(log_rows).to_csv(log_csv, index=False)

    print(f"[DONE] Connectivity matrix: {output_csv} | shape={df.shape}")
    print(f"[DONE] Processing log: {log_csv}")


if __name__ == "__main__":
    main()
