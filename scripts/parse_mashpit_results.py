#!/usr/bin/env python3
"""Normalize supported Mashpit candidate files into a stable result schema."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from common import WorkflowError, write_json


CLUSTER_FIELDS = ("PDS_acc", "cluster", "cluster_id", "snp_cluster", "pd_cluster")
SCORE_FIELDS = (
    "best_similarity_score", "similarity_score", "score", "similarity", "jaccard"
)


def first_value(row: dict[str, Any], names: tuple[str, ...]) -> Any:
    for name in names:
        if row.get(name) not in (None, ""):
            return row[name]
    return None


def normalize_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for row in rows:
        cluster = first_value(row, CLUSTER_FIELDS)
        score_value = first_value(row, SCORE_FIELDS)
        if cluster is None or score_value is None:
            continue
        try:
            score = float(score_value)
        except (TypeError, ValueError):
            continue
        candidates.append({
            "cluster": str(cluster),
            "score": score,
            "near_top": str(row.get("near_top", "")).lower() in {"true", "1", "yes"},
            "metadata": row,
        })
    candidates.sort(key=lambda item: item["score"], reverse=True)
    return candidates


def load_candidates(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() == ".json":
        value = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(value, dict):
            value = value.get("candidates", value.get("results", []))
        if not isinstance(value, list):
            raise WorkflowError("Mashpit JSON must contain a candidate list.")
        return normalize_rows([row for row in value if isinstance(row, dict)])
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(encoding="utf-8", newline="") as handle:
        return normalize_rows(list(csv.DictReader(handle, delimiter=delimiter)))


def locate_candidate_file(output_dir: Path) -> Path:
    matches = list(output_dir.glob("*_cluster_candidates.csv"))
    if len(matches) != 1:
        raise WorkflowError("Expected exactly one Mashpit *_cluster_candidates.csv file.")
    return matches[0]


def locate_representative_file(output_dir: Path) -> Path:
    matches = list(output_dir.glob("*_representative_matches.csv"))
    if len(matches) != 1:
        raise WorkflowError("Expected exactly one Mashpit *_representative_matches.csv file.")
    return matches[0]


def locate_tree_image(output_dir: Path) -> Path | None:
    # Unlike the two locators above, a missing tree is not an error: Mashpit
    # itself skips tree generation when the top hit is below --threshold or
    # fewer than two candidates qualify (see references/mashpit-interpretation.md).
    matches = list(output_dir.glob("*_tree.png"))
    return matches[0] if len(matches) == 1 else None


def best_representative_score(output_dir: Path) -> float:
    path = locate_representative_file(output_dir)
    with path.open(encoding="utf-8", newline="") as handle:
        scores = []
        for row in csv.DictReader(handle):
            try:
                scores.append(float(row["similarity_score"]))
            except (KeyError, TypeError, ValueError):
                continue
    if not scores:
        raise WorkflowError("Mashpit representative table has no usable similarity_score values.")
    return max(scores)


def interpret(candidates: list[dict[str, Any]], threshold: float | None = None) -> dict:
    if not candidates:
        return {
            "status": "NO_MATCH",
            "best_candidate": None,
            "alternative_candidates": [],
            "ambiguous": False,
            "below_threshold": False,
            "screening_result": "No candidate",
            "warnings": ["Mashpit returned no candidate; there was no usable sketch overlap with the query."],
        }
    best = candidates[0]
    tied = [item for item in candidates[1:] if item.get("near_top")]
    ambiguous = bool(tied)
    # --number/--tie-tolerance-hashes return and group the top hits regardless of how low
    # their absolute score is; --threshold only gates local tree construction (see
    # references/mashpit-interpretation.md). A "near-top"/single-cluster result can
    # therefore still be noise-level similarity, so that has to be checked separately.
    below_threshold = threshold is not None and best["score"] < threshold
    warnings = []
    if ambiguous:
        warnings.append("Mashpit marks multiple SNP clusters as near-top within sketch resolution.")
    if below_threshold:
        warnings.append(
            f"Mashpit's top hit (score {best['score']:.3f}) is below its own query threshold "
            f"({threshold}); treat this as noise-level sketch similarity, not a real candidate."
        )
    return {
        "status": "AMBIGUOUS" if ambiguous else "MATCH",
        "best_candidate": {"cluster": best["cluster"], "score": best["score"]},
        "alternative_candidates": [
            {"cluster": item["cluster"], "score": item["score"]} for item in candidates[1:6]
        ],
        "ambiguous": ambiguous,
        "below_threshold": below_threshold,
        "screening_result": "Ambiguous top candidates" if ambiguous else "Top-ranked candidate",
        "warnings": warnings,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    parser.add_argument("--mashpit-output-dir")
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        source = Path(args.input) if args.input else locate_candidate_file(Path(args.mashpit_output_dir))
        result = interpret(load_candidates(source), threshold=args.threshold)
        result["source_file"] = str(source)
    except (WorkflowError, OSError, ValueError, json.JSONDecodeError) as error:
        result = {"status": "ERROR", "error": str(error), "best_candidate": None, "ambiguous": True}
    write_json(Path(args.output), result)
    return 0 if result["status"] != "ERROR" else 2


if __name__ == "__main__":
    raise SystemExit(main())
