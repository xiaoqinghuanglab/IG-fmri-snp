#!/usr/bin/env python3
"""Compute MSDL (39 ROI) functional connectivity features from preprocessed rs-fMRI.

What it does (high level):
  - loads each subject's 4D NIfTI
  - extracts ROI time series using the MSDL atlas
  - optionally aligns all subjects to a fixed number of time points (truncate/pad)
  - fits GraphicalLassoCV and uses the precision matrix as sparse connectivity
  - writes a subject-by-edge CSV with 741 edges (39*38/2)

Inputs
  --input-dir: directory of 4D NIfTI files (*.nii.gz)
  --msdl-maps: optional local atlas maps file; if omitted, will try nilearn's fetch_atlas_msdl()

Output
  --output-csv: CSV with:
    - Subject_ID column
    - Connectivity_<ROI1>_<ROI2> columns

Example
  python scripts/fmri/compute_connectivity_msdl.py \
    --input-dir /path/processed_nifti \
    --output-csv outputs/connectivity/Connectivity_matrix_all_subjects_region_pairs.csv \
    --tr 3 \
    --fixed-length median
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

import nibabel as nib
from sklearn.covariance import GraphicalLassoCV

try:
    from nilearn.datasets import fetch_atlas_msdl
    from nilearn.maskers import NiftiMapsMasker
except Exception as e:
    raise SystemExit(
        "nilearn is required for this script. Install dependencies via: pip install -r requirements.txt"
    ) from e


def parse_args():
    p = argparse.ArgumentParser(description="Compute MSDL sparse connectivity features (GraphicalLassoCV).")
    p.add_argument("--input-dir", type=Path, required=True, help="Folder containing preprocessed 4D NIfTI files.")
    p.add_argument("--pattern", type=str, default="*_smooth.nii.gz", help="Glob pattern for input NIfTI files.")
    p.add_argument("--output-csv", type=Path, required=True, help="Output CSV path (subject-by-edge matrix).")
    p.add_argument("--tr", type=float, default=3.0, help="TR in seconds (used by band-pass filtering).")
    p.add_argument("--high-pass", type=float, default=0.01, help="High-pass cutoff (Hz).")
    p.add_argument("--low-pass", type=float, default=0.1, help="Low-pass cutoff (Hz).")
    p.add_argument("--detrend", action="store_true", help="Detrend ROI time series.")
    p.add_argument("--no-standardize", action="store_true", help="Disable standardization (z-scoring).")
    p.add_argument("--msdl-maps", type=Path, default=None, help="Optional local path to the MSDL atlas maps NIfTI.")
    p.add_argument("--fixed-length", type=str, default="none",
                   choices=["none", "median", "min"],
                   help="Time-series length alignment: none|min|median (truncate/pad).")
    p.add_argument("--id-regex", type=str, default=r"(\d{3}_S_\d{4})",
                   help="Regex with one capture group to extract Subject_ID from filename.")
    p.add_argument("--abs-edges", action="store_true", help="Use absolute precision values for edges (default in paper).")
    p.add_argument("--random-state", type=int, default=0, help="Random state for GraphicalLassoCV.")
    return p.parse_args()


def subject_id_from_path(p: Path, id_re: re.Pattern) -> str:
    m = id_re.search(p.name)
    if m:
        return m.group(1)
    return p.stem.replace(".nii", "")


def infer_T(img_path: Path) -> int:
    img = nib.load(str(img_path))
    shape = img.shape
    if len(shape) < 4:
        raise ValueError(f"Expected 4D image, got shape={shape} for {img_path}")
    return int(shape[3])


def align_length(ts: np.ndarray, L: int) -> np.ndarray:
    """Truncate or zero-pad time series (T x R) to length L."""
    T, R = ts.shape
    if T == L:
        return ts
    if T > L:
        return ts[:L, :]
    out = np.zeros((L, R), dtype=ts.dtype)
    out[:T, :] = ts
    return out


def main():
    args = parse_args()
    in_dir = args.input_dir.expanduser().resolve()
    files = sorted(in_dir.glob(args.pattern))
    if not files:
        raise SystemExit(f"No files found in {in_dir} with pattern '{args.pattern}'")

    id_re = re.compile(args.id_regex)

    # Load atlas
    if args.msdl_maps is None:
        msdl = fetch_atlas_msdl()
        maps_img = msdl["maps"]
        labels = [str(x) for x in msdl["labels"]]
    else:
        maps_img = str(args.msdl_maps.expanduser().resolve())
        # If you provide custom maps, also provide labels separately if needed.
        # We'll fall back to ROI_0..ROI_38
        labels = [f"ROI_{i}" for i in range(39)]

    masker = NiftiMapsMasker(
        maps_img=maps_img,
        detrend=args.detrend,
        standardize=not args.no_standardize,
        high_pass=args.high_pass,
        low_pass=args.low_pass,
        t_r=args.tr,
    )

    # Decide fixed length
    target_L = None
    if args.fixed_length != "none":
        lengths = [infer_T(f) for f in files]
        if args.fixed_length == "min":
            target_L = int(np.min(lengths))
        elif args.fixed_length == "median":
            target_L = int(np.median(lengths))

    rows = []
    for f in files:
        sid = subject_id_from_path(f, id_re)

        ts = masker.fit_transform(str(f))  # (T, 39)
        if target_L is not None:
            ts = align_length(ts, target_L)

        # Graphical lasso -> precision matrix
        gl = GraphicalLassoCV()
        gl.fit(ts)
        precision = gl.precision_.copy()

        # build edges (upper triangle, exclude diagonal)
        n = precision.shape[0]
        feats = {}
        for i in range(n):
            for j in range(i + 1, n):
                v = precision[i, j]
                if args.abs_edges:
                    v = float(abs(v))
                else:
                    v = float(v)
                feats[f"Connectivity_{labels[i]}_{labels[j]}"] = v

        feats["Subject_ID"] = sid
        rows.append(feats)

    df = pd.DataFrame(rows)
    # stable column order: Subject_ID first, then edges sorted
    edge_cols = sorted([c for c in df.columns if c != "Subject_ID"])
    df = df[["Subject_ID"] + edge_cols]

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    print(f"[DONE] Wrote: {args.output_csv} | shape={df.shape}")


if __name__ == "__main__":
    main()
