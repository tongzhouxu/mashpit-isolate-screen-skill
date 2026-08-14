#!/usr/bin/env python3
"""Select query-relevant representative genomes for ska2 SNP resolution."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from common import WorkflowError, load_json, write_json
from parse_mashpit_results import (
    CLUSTER_FIELDS,
    SCORE_FIELDS,
    first_value,
    load_candidates,
    locate_candidate_file,
    locate_representative_file,
)


ACCESSION_FIELDS = ("asm_acc", "assembly_accession", "assembly_acc")
BIOSAMPLE_FIELDS = ("biosample_acc", "biosample_accession")


def relevant_clusters(mashpit_output_dir: Path) -> list[str]:
    candidates = load_candidates(locate_candidate_file(mashpit_output_dir))
    if not candidates:
        return []
    clusters = [candidates[0]["cluster"]]
    clusters.extend(item["cluster"] for item in candidates[1:] if item.get("near_top"))
    return clusters


def load_representatives(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    representatives = []
    for row in rows:
        cluster = first_value(row, CLUSTER_FIELDS)
        accession = first_value(row, ACCESSION_FIELDS)
        score_value = first_value(row, SCORE_FIELDS)
        if cluster is None or accession is None or score_value is None:
            continue
        try:
            score = float(score_value)
        except (TypeError, ValueError):
            continue
        representatives.append({
            "cluster": str(cluster),
            "accession": str(accession),
            "biosample": first_value(row, BIOSAMPLE_FIELDS),
            "score": score,
        })
    return representatives


def select_targets(mashpit_output_dir: Path, policy: dict[str, Any]) -> dict[str, Any]:
    clusters = relevant_clusters(mashpit_output_dir)
    if not clusters:
        return {"status": "SKIPPED", "reason": "No Mashpit candidate to resolve.", "targets": []}
    representatives = load_representatives(locate_representative_file(mashpit_output_dir))
    by_cluster: dict[str, list[dict[str, Any]]] = {cluster: [] for cluster in clusters}
    for row in representatives:
        if row["cluster"] in by_cluster:
            by_cluster[row["cluster"]].append(row)
    max_per_cluster = policy["max_representatives_per_cluster"]
    for cluster in by_cluster:
        by_cluster[cluster].sort(key=lambda item: item["score"], reverse=True)
        by_cluster[cluster] = by_cluster[cluster][:max_per_cluster]
    if not any(by_cluster.values()):
        return {
            "status": "SKIPPED",
            "reason": "No representative genomes with usable accessions in the relevant clusters.",
            "targets": [],
        }
    # Round-robin across relevant clusters so each one keeps at least one
    # representative even after the total-genome cap is applied.
    ordered: list[dict[str, Any]] = []
    max_len = max((len(rows) for rows in by_cluster.values()), default=0)
    max_total = policy["max_total_genomes"]
    round_index = 0
    while round_index < max_len and len(ordered) < max_total:
        for cluster in clusters:
            rows = by_cluster[cluster]
            if round_index < len(rows):
                ordered.append(rows[round_index])
        round_index += 1
    return {
        "status": "SELECTED",
        "relevant_clusters": clusters,
        "targets": ordered[:max_total],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mashpit-output-dir", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    policy = load_json(Path(args.policy))
    try:
        result = select_targets(Path(args.mashpit_output_dir), policy)
    except (WorkflowError, OSError, ValueError) as error:
        result = {"status": "ERROR", "error": str(error), "targets": []}
    write_json(Path(args.output), result)
    return 0 if result["status"] != "ERROR" else 2


if __name__ == "__main__":
    raise SystemExit(main())
