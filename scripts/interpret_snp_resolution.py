#!/usr/bin/env python3
"""Rank ska2 pairwise SNP distances between the query and candidate representatives."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from common import write_json


QUERY_SAMPLE_ID = "QUERY"


def query_distances(distance_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    results = []
    for row in distance_rows:
        if row["sample1"] == QUERY_SAMPLE_ID:
            other = row["sample2"]
        elif row["sample2"] == QUERY_SAMPLE_ID:
            other = row["sample1"]
        else:
            continue
        results.append({
            "sample": other,
            "snp_distance": row["snp_distance"],
            "mismatch_proportion": row["mismatch_proportion"],
            "match_count": row["match_count"],
            "mismatch_count": row["mismatch_count"],
        })
    results.sort(key=lambda item: item["snp_distance"])
    return results


def interpret(
    distance_rows: list[dict[str, Any]],
    targets: list[dict[str, Any]],
    mash_best_cluster: str | None,
) -> dict[str, Any]:
    ranked = query_distances(distance_rows)
    if not ranked:
        return {
            "status": "INSUFFICIENT_DATA",
            "ranked": [],
            "warnings": ["No query-to-representative SNP distances were computed."],
        }
    accession_to_cluster = {item["accession"]: item["cluster"] for item in targets}
    for item in ranked:
        item["cluster"] = accession_to_cluster.get(item["sample"])
    nearest = ranked[0]
    agrees_with_mash = mash_best_cluster is not None and nearest["cluster"] == mash_best_cluster
    warnings = []
    if mash_best_cluster is not None and not agrees_with_mash:
        warnings.append(
            "The SNP-nearest genome belongs to a different cluster "
            f"({nearest['cluster']}) than Mashpit's top candidate ({mash_best_cluster})."
        )
    return {
        "status": "RESOLVED",
        "nearest_sample": nearest["sample"],
        "nearest_cluster": nearest["cluster"],
        "nearest_snp_distance": nearest["snp_distance"],
        "agrees_with_mash_top_candidate": agrees_with_mash,
        "ranked": ranked,
        "warnings": warnings,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--distances-json", required=True)
    parser.add_argument("--targets-json", required=True)
    parser.add_argument("--mash-best-cluster")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    distances = json.loads(Path(args.distances_json).read_text(encoding="utf-8"))
    targets = json.loads(Path(args.targets_json).read_text(encoding="utf-8"))
    result = interpret(distances, targets, args.mash_best_cluster)
    write_json(Path(args.output), result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
