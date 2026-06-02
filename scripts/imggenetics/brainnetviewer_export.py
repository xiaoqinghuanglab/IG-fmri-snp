#!/usr/bin/env python3
"""Export BrainNetViewer node and edge files from pairwise SNP-connectivity results."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


MSDL_LABELS = [
    "Auditory_L",
    "Auditory_R",
    "Striate",
    "DMN_L",
    "DMN_Med",
    "DMN_Front",
    "DMN_R",
    "Occ_post",
    "Motor",
    "DLPFC_R",
    "Front_pol_R",
    "Par_R",
    "Post_Temp_R",
    "Basal",
    "Par_L",
    "DLPFC_L",
    "Front_pol_L",
    "IPS_L",
    "IPS_R",
    "LOC_L",
    "Vis",
    "LOC_R",
    "D_ACC",
    "V_ACC",
    "A_Ins_R",
    "STS_L",
    "STS_R",
    "TPJ_L",
    "Broca",
    "Sup_Front_S",
    "TPJ_R",
    "Pars_Op_R",
    "Cereb",
    "Dors_PCC",
    "Ins_L",
    "Cing",
    "Ins_R",
    "Ant_IPS_L",
    "Ant_IPS_R",
]

MSDL_COORDS = np.array(
    [
        [-52.5, -22.5, 7.5],
        [52.5, -22.5, 7.5],
        [-2.5, -82.5, 17.5],
        [-57.5, -42.5, 27.5],
        [-2.5, -57.5, 37.5],
        [-2.5, 47.5, 17.5],
        [57.5, -47.5, 27.5],
        [-2.5, -82.5, 47.5],
        [-2.5, -22.5, 57.5],
        [42.5, 37.5, 37.5],
        [27.5, 62.5, 7.5],
        [47.5, -67.5, 47.5],
        [52.5, -42.5, 7.5],
        [-2.5, -12.5, -7.5],
        [-47.5, -67.5, 47.5],
        [-42.5, 37.5, 37.5],
        [-27.5, 62.5, 7.5],
        [-42.5, -67.5, 37.5],
        [42.5, -67.5, 37.5],
        [-42.5, -77.5, 7.5],
        [-2.5, -92.5, -7.5],
        [42.5, -77.5, 7.5],
        [-2.5, 27.5, 27.5],
        [-2.5, 12.5, 7.5],
        [37.5, 12.5, 7.5],
        [-52.5, -22.5, -7.5],
        [52.5, -22.5, -7.5],
        [-52.5, -57.5, 7.5],
        [-42.5, 17.5, 27.5],
        [-2.5, 57.5, 37.5],
        [52.5, -57.5, 7.5],
        [42.5, 17.5, 27.5],
        [-2.5, -67.5, -27.5],
        [-2.5, -47.5, 47.5],
        [-37.5, 12.5, 7.5],
        [-2.5, 7.5, 27.5],
        [37.5, 12.5, 7.5],
        [-37.5, -67.5, 27.5],
        [37.5, -67.5, 27.5],
    ]
)

LABEL_TO_IDX = {label: idx for idx, label in enumerate(MSDL_LABELS)}


def parse_args():
    parser = argparse.ArgumentParser(description="Export BrainNetViewer edge lists from pairwise imggenetics results.")
    parser.add_argument("--results-csv", type=Path, required=True, help="Per-comparison results CSV.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument(
        "--edge-col",
        type=str,
        default="edge_id",
        help="Column containing connectivity edge identifiers (default: edge_id).",
    )
    parser.add_argument(
        "--weight-col",
        type=str,
        default="beta_interaction",
        help="Column to use as edge weights (default: beta_interaction).",
    )
    parser.add_argument(
        "--rank-col",
        type=str,
        default="fdr_interaction",
        help="Column used to keep the lowest-ranked row per edge (default: fdr_interaction).",
    )
    parser.add_argument(
        "--label",
        type=str,
        default=None,
        help="Optional label for the edge filename. Defaults to the input stem.",
    )
    return parser.parse_args()


def parse_edge_name(name: str):
    text = str(name)
    if not text.startswith("Connectivity_"):
        return None
    parts = text.split("_")
    if len(parts) < 3:
        return None
    roi1 = parts[1]
    roi2 = "_".join(parts[2:])
    return roi1, roi2


def write_node_file(path: Path):
    lines = []
    for (x, y, z), label in zip(MSDL_COORDS, MSDL_LABELS):
        lines.append(f"{x}\t{y}\t{z}\t1\t{label}")
    path.write_text("\n".join(lines) + "\n")


def main():
    args = parse_args()
    df = pd.read_csv(args.results_csv)

    for required_col in [args.edge_col, args.weight_col]:
        if required_col not in df.columns:
            raise SystemExit(f"Missing required column '{required_col}' in: {args.results_csv}")

    if args.rank_col in df.columns:
        df["_rank"] = pd.to_numeric(df[args.rank_col], errors="coerce")
        df = df.sort_values("_rank", ascending=True, na_position="last")

    df = df.drop_duplicates(subset=[args.edge_col], keep="first").copy()
    edges = []

    for _, row in df.iterrows():
        parsed = parse_edge_name(row[args.edge_col])
        if not parsed:
            continue

        roi1, roi2 = parsed
        if roi1 not in LABEL_TO_IDX or roi2 not in LABEL_TO_IDX:
            continue

        i = LABEL_TO_IDX[roi1] + 1
        j = LABEL_TO_IDX[roi2] + 1
        weight = float(row[args.weight_col])
        edges.append((i, j, weight))

    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    label = args.label or args.results_csv.stem
    node_path = out_dir / "msdl_nodes.node"
    edge_path = out_dir / f"edges_{label}.edge"
    write_node_file(node_path)
    pd.DataFrame(edges, columns=["i", "j", "weight"]).to_csv(edge_path, sep="\t", header=False, index=False)

    print(f"[DONE] nodes: {node_path}")
    print(f"[DONE] edges: {edge_path} | n_edges={len(edges)}")


if __name__ == "__main__":
    main()
