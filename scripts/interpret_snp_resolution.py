#!/usr/bin/env python3
"""Interpret ska2 pairwise SNP distances: rank, group by cluster, and build a tree.

Reports the full pairwise matrix ska2 already computes (not just query-vs-
representative), a per-cluster distance summary, a Neighbor-Joining tree over
every compared genome, and a structural confidence comparison (nearest
cluster's closest genome vs. the next-nearest cluster's closest genome).
No strength label ("high/medium/low confidence") is invented: only the raw
distances and their ratio are reported, consistent with
references/mashpit-interpretation.md's policy on the Mash score itself.
"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path
from typing import Any

from build_snp_tree import neighbor_joining
from common import WorkflowError, write_json


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


def _display_label(sample: str, accession_to_cluster: dict[str, str]) -> str:
    cluster = accession_to_cluster.get(sample)
    return sample if cluster is None else f"{sample}_{cluster}"


def _relabel_rows(
    distance_rows: list[dict[str, Any]], accession_to_cluster: dict[str, str]
) -> list[dict[str, Any]]:
    relabeled = []
    for row in distance_rows:
        relabeled.append({
            **row,
            "sample1": _display_label(row["sample1"], accession_to_cluster),
            "sample2": _display_label(row["sample2"], accession_to_cluster),
        })
    return relabeled


def cluster_summary(ranked: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_cluster: dict[str, list[float]] = {}
    for item in ranked:
        by_cluster.setdefault(item["cluster"], []).append(item["snp_distance"])
    summary = [
        {
            "cluster": cluster,
            "genomes_compared": len(distances),
            "min_snp_distance": min(distances),
            "median_snp_distance": statistics.median(distances),
            "max_snp_distance": max(distances),
        }
        for cluster, distances in by_cluster.items()
    ]
    summary.sort(key=lambda item: item["min_snp_distance"])
    return summary


def build_confidence(summary: list[dict[str, Any]]) -> dict[str, Any]:
    nearest = summary[0]
    result: dict[str, Any] = {
        "nearest_cluster": nearest["cluster"],
        "nearest_cluster_min_snp_distance": nearest["min_snp_distance"],
        "next_cluster": None,
        "next_cluster_min_snp_distance": None,
        "separation_ratio": None,
    }
    if len(summary) == 1:
        result["statement"] = (
            f"Only one cluster ({nearest['cluster']}) had representatives close enough for Mashpit to "
            f"flag as plausible; ska2 independently confirms it at {nearest['min_snp_distance']:.2f} SNPs. "
            "There is no second candidate to compare against at SNP resolution."
        )
        return result
    runner_up = summary[1]
    result["next_cluster"] = runner_up["cluster"]
    result["next_cluster_min_snp_distance"] = runner_up["min_snp_distance"]
    if nearest["min_snp_distance"] == 0:
        result["statement"] = (
            f"The query is a 0-SNP match within cluster {nearest['cluster']}, versus "
            f"{runner_up['min_snp_distance']:.2f} SNPs to the nearest genome in the next-closest cluster "
            f"({runner_up['cluster']}) — an exact match, not merely a nearby one."
        )
        return result
    ratio = runner_up["min_snp_distance"] / nearest["min_snp_distance"]
    result["separation_ratio"] = ratio
    result["statement"] = (
        f"The nearest genome in cluster {nearest['cluster']} is {nearest['min_snp_distance']:.2f} SNPs away, "
        f"versus {runner_up['min_snp_distance']:.2f} SNPs to the nearest genome in the next-closest cluster "
        f"({runner_up['cluster']}) — about {ratio:.1f}x closer. This is a raw distance ratio, not a "
        "statistical confidence score; no validated SNP cutoff for cluster assignment has been supplied."
    )
    return result


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
    summary = cluster_summary(ranked)
    try:
        newick = neighbor_joining(_relabel_rows(distance_rows, accession_to_cluster))
        tree_status = "PASS"
    except WorkflowError as error:
        newick = None
        tree_status = str(error)
    return {
        "status": "RESOLVED",
        "nearest_sample": nearest["sample"],
        "nearest_cluster": nearest["cluster"],
        "nearest_snp_distance": nearest["snp_distance"],
        "agrees_with_mash_top_candidate": agrees_with_mash,
        "ranked": ranked,
        "cluster_summary": summary,
        "confidence": build_confidence(summary),
        "pairwise_distances": distance_rows,
        "newick_tree": newick,
        "tree_status": tree_status,
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
