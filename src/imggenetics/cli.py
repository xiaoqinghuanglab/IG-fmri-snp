"""Command-line entry point for the paper-aligned pairwise analysis workflow."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

from .pairwise import run_pairwise_analysis


def first_existing_path(candidates: Iterable[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def resolve_input_path(explicit: Path | None, candidates: list[Path], label: str) -> Path:
    if explicit is not None:
        return explicit

    resolved = first_existing_path(candidates)
    if resolved is None:
        tried = "\n".join(f"  - {candidate}" for candidate in candidates)
        raise FileNotFoundError(f"Could not resolve {label}. Tried:\n{tried}")
    return resolved


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="imggenetics",
        description="Pairwise SNP-connectivity interaction analysis for the IG-fmri-snp paper workflow.",
    )

    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data"),
        help="Folder containing the standard pipeline CSV inputs.",
    )
    parser.add_argument("--genetic-file", type=Path, default=None, help="Union genotype matrix CSV.")
    parser.add_argument("--connectivity-file", type=Path, default=None, help="MSDL subject-by-edge connectivity CSV.")
    parser.add_argument("--covariates-file", type=Path, default=None, help="Covariate CSV with subtype labels.")
    parser.add_argument("--snp-metadata-file", type=Path, default=None, help="Combined SNP metadata CSV.")
    parser.add_argument(
        "--candidate-edges-file",
        type=Path,
        default=None,
        help="model2 candidate resilience edge CSV from scripts/fmriphenotype/get_resilience_edge.py",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs") / "imggenetics",
        help="Folder where per-comparison results will be written.",
    )
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    data_dir = args.data_dir
    phenotype_output_dir = Path("outputs") / "fmriphenotype"
    connectivity_path = resolve_input_path(
        args.connectivity_file,
        [
            data_dir / "msdl_all_subjects_connectivity_edges.csv",
            data_dir / "Connectivity_matrix_all_subjects_region_pairs.csv",
        ],
        "connectivity file",
    )
    covariates_path = resolve_input_path(
        args.covariates_file,
        [
            data_dir / "covariate_file.csv",
            data_dir / "subject_covariates.csv",
        ],
        "covariates file",
    )
    genetic_path = resolve_input_path(
        args.genetic_file,
        [data_dir / "genotype_matrix.csv"],
        "genotype matrix",
    )
    snp_metadata_path = resolve_input_path(
        args.snp_metadata_file,
        [data_dir / "snp_metadata.csv"],
        "SNP metadata file",
    )
    candidate_edge_candidates = [
        data_dir / "model2_candidate_resilience_edges.csv",
        data_dir / "Model2_MSDL_Phenotype" / "model2_candidate_resilience_edges.csv",
        phenotype_output_dir / "model2_candidate_resilience_edges.csv",
    ]
    try:
        candidate_edges_path = resolve_input_path(
            args.candidate_edges_file,
            candidate_edge_candidates,
            "candidate resilience edge file",
        )
    except FileNotFoundError as exc:
        tried = "\n".join(f"  - {candidate}" for candidate in candidate_edge_candidates)
        raise FileNotFoundError(
            "Could not resolve the candidate resilience edge file.\n"
            "Run `bash run.sh python scripts/fmriphenotype/get_resilience_edge.py` first,\n"
            "or pass an explicit --candidate-edges-file path.\n"
            f"Tried:\n{tried}"
        ) from exc

    run_pairwise_analysis(
        connectivity_path=connectivity_path,
        covariates_path=covariates_path,
        genetic_path=genetic_path,
        snp_metadata_path=snp_metadata_path,
        candidate_edges_path=candidate_edges_path,
        output_dir=args.output_dir,
    )
    print(f"Done. Outputs written to: {args.output_dir.resolve()}")
    return 0
