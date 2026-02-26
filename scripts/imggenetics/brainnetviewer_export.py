#!/usr/bin/env python3
"""Create BrainNet Viewer node/edge files from significant SNP-edge results.

This utility focuses on connectivity visualization:
  - aggregates significant SNP-edge pairs to unique edges (keeping the smallest FDR if available)
  - exports a BrainNet edge list using MSDL ROI coordinates (39 ROI)

Inputs:
  --significant-csv : outputs/significant_result.csv from imggenetics
  --contrast        : which P_FDR / Coeff columns to use (e.g., AsymAD_vs_Control)
Outputs:
  --out-dir:
    - msdl_nodes.node  (ROI coordinates + labels)
    - edges_<contrast>.edge (edge list)

Notes:
  - BrainNet Viewer has multiple edge formats; this script writes a simple 3-column edge list:
        i  j  weight
    where i/j are 1-based ROI indices in the node file.
  - Connectivity_Name must encode ROI labels in a consistent pattern:
        Connectivity_<ROI1>_<ROI2>
    (this is the default output from compute_connectivity_msdl.py)

Example:
  python scripts/imggenetics/brainnetviewer_export.py \
    --significant-csv outputs/significant_result.csv \
    --contrast AsymAD_vs_Control \
    --out-dir outputs/brainnet
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# MSDL atlas (39 ROI) coordinates + labels (as used in the notebook).
# If you prefer an external node file, you can replace this section.
MSDL_LABELS = ['Auditory_L', 'Auditory_R', 'Striate', 'DMN_L', 'DMN_Med', 'DMN_Front', 'DMN_R', 'Occ_post',
              'Motor', 'DLPFC_R', 'Front_pol_R', 'Par_R', 'Post_Temp_R', 'Basal', 'Par_L', 'DLPFC_L',
              'Front_pol_L', 'IPS_L', 'IPS_R', 'LOC_L', 'Vis', 'LOC_R', 'D_ACC', 'V_ACC', 'A_Ins_R', 'STS_L',
              'STS_R', 'TPJ_L', 'Broca', 'Sup_Front_S', 'TPJ_R', 'Pars_Op_R', 'Cereb', 'Dors_PCC', 'Ins_L',
              'Cing', 'Ins_R', 'Ant_IPS_L', 'Ant_IPS_R']

MSDL_COORDS = np.array([
    [-52.5, -22.5, 7.5], [52.5, -22.5, 7.5], [-2.5, -82.5, 17.5], [-57.5, -42.5, 27.5], [-2.5, -57.5, 37.5],
    [-2.5, 47.5, 17.5], [57.5, -47.5, 27.5], [-2.5, -82.5, 47.5], [-2.5, -22.5, 57.5], [42.5, 37.5, 37.5],
    [27.5, 62.5, 7.5], [47.5, -67.5, 47.5], [52.5, -42.5, 7.5], [-2.5, -12.5, -7.5], [-47.5, -67.5, 47.5],
    [-42.5, 37.5, 37.5], [-27.5, 62.5, 7.5], [-42.5, -67.5, 37.5], [42.5, -67.5, 37.5], [-42.5, -77.5, 7.5],
    [-2.5, -92.5, -7.5], [42.5, -77.5, 7.5], [-2.5, 27.5, 27.5], [-2.5, 12.5, 7.5], [37.5, 12.5, 7.5],
    [-52.5, -22.5, -7.5], [52.5, -22.5, -7.5], [-52.5, -57.5, 7.5], [-42.5, 17.5, 27.5], [-2.5, 57.5, 37.5],
    [52.5, -57.5, 7.5], [42.5, 17.5, 27.5], [-2.5, -67.5, -27.5], [-2.5, -47.5, 47.5], [-37.5, 12.5, 7.5],
    [-2.5, 7.5, 27.5], [37.5, 12.5, 7.5], [-37.5, -67.5, 27.5], [37.5, -67.5, 27.5],
])

LABEL_TO_IDX = {lab: i for i, lab in enumerate(MSDL_LABELS)}  # 0-based


def parse_args():
    p = argparse.ArgumentParser(description="Export BrainNet Viewer files from significant SNP-edge results.")
    p.add_argument("--significant-csv", type=Path, required=True, help="significant_result.csv from imggenetics.")
    p.add_argument("--contrast", type=str, required=True, help="Contrast label used in column names (e.g., AsymAD_vs_Control).")
    p.add_argument("--out-dir", type=Path, required=True, help="Output folder.")
    p.add_argument("--fdr-col", type=str, default=None, help="Override FDR column name (otherwise inferred).")
    p.add_argument("--coeff-col", type=str, default=None, help="Override coefficient column name (otherwise inferred).")
    return p.parse_args()


def parse_edge_name(name: str):
    # Expect Connectivity_<ROI1>_<ROI2>
    s = str(name)
    if not s.startswith("Connectivity_"):
        return None
    parts = s.split("_")
    if len(parts) < 3:
        return None
    roi1 = parts[1]
    roi2 = "_".join(parts[2:])
    return roi1, roi2


def write_node_file(path: Path):
    # BrainNet node format commonly: x y z size label
    # Here we use size=1 for all nodes
    lines = []
    for (x, y, z), lab in zip(MSDL_COORDS, MSDL_LABELS):
        lines.append(f"{x}\t{y}\t{z}\t1\t{lab}")
    path.write_text("\n".join(lines) + "\n")


def main():
    args = parse_args()
    df = pd.read_csv(args.significant_csv)

    fdr_col = args.fdr_col or f"P_FDR_{args.contrast}"
    coeff_col = args.coeff_col or f"Coeff_{args.contrast}"

    if fdr_col not in df.columns:
        # Some pipelines use P_FDR<...> without underscore after FDR
        alt = f"P_FDR{args.contrast}"
        if alt in df.columns:
            fdr_col = alt
    if coeff_col not in df.columns:
        raise SystemExit(f"Coefficient column not found: {coeff_col}")

    # Keep rows that have a coeff for this contrast
    sub = df.dropna(subset=[coeff_col]).copy()
    if sub.empty:
        raise SystemExit("No rows with coefficients for the requested contrast.")

    # Collapse to unique edges (pick smallest FDR if available, else first)
    if fdr_col in sub.columns:
        sub["_rank"] = pd.to_numeric(sub[fdr_col], errors="coerce")
    else:
        sub["_rank"] = np.nan

    sub = sub.sort_values("_rank", ascending=True, na_position="last")
    sub = sub.drop_duplicates(subset=["Connectivity_Name"], keep="first")

    # Build edge list
    edges = []
    for _, r in sub.iterrows():
        parsed = parse_edge_name(r["Connectivity_Name"])
        if not parsed:
            continue
        roi1, roi2 = parsed
        if roi1 not in LABEL_TO_IDX or roi2 not in LABEL_TO_IDX:
            continue
        i = LABEL_TO_IDX[roi1] + 1
        j = LABEL_TO_IDX[roi2] + 1
        w = float(r[coeff_col])
        edges.append((i, j, w))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    node_path = args.out_dir / "msdl_nodes.node"
    edge_path = args.out_dir / f"edges_{args.contrast}.edge"
    write_node_file(node_path)

    edf = pd.DataFrame(edges, columns=["i", "j", "weight"])
    edf.to_csv(edge_path, sep="\t", header=False, index=False)

    print(f"[DONE] nodes: {node_path}")
    print(f"[DONE] edges: {edge_path} | n_edges={len(edf)}")


if __name__ == "__main__":
    main()
