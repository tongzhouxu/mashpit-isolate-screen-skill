#!/usr/bin/env python3
"""Build a Neighbor-Joining tree (Newick) from a pairwise ska2 SNP distance table.

Standard Saitou-Nei neighbor joining, implemented directly rather than adding a
new pinned tree-building tool: the input is already a small, fully-resolved
pairwise distance matrix (<= max_total_genomes taxa), which NJ handles exactly
and cheaply. Verified against a hand-computed additive 4-taxon case in
tests/test_workflow.py, where NJ is guaranteed to recover the true tree exactly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from common import WorkflowError, write_json


def _matrix_from_rows(rows: list[dict[str, Any]]) -> tuple[list[str], dict[tuple[str, str], float]]:
    labels = sorted({row["sample1"] for row in rows} | {row["sample2"] for row in rows})
    distances: dict[tuple[str, str], float] = {}
    for row in rows:
        a, b, d = row["sample1"], row["sample2"], float(row["snp_distance"])
        distances[(a, b)] = d
        distances[(b, a)] = d
    return labels, distances


def neighbor_joining(rows: list[dict[str, Any]]) -> str:
    labels, d = _matrix_from_rows(rows)
    if len(labels) < 2:
        raise WorkflowError("Neighbor joining requires at least two samples.")
    if len(labels) == 2:
        a, b = labels
        dab = d[(a, b)]
        return f"({a}:{dab / 2:.4f},{b}:{dab / 2:.4f});"

    nodes: dict[str, str] = {label: label for label in labels}
    active = list(labels)
    counter = 0

    while len(active) > 2:
        n = len(active)
        total = {i: sum(d[(i, k)] for k in active if k != i) for i in active}
        best = None
        best_q = None
        for a_idx in range(n):
            for b_idx in range(a_idx + 1, n):
                i, j = active[a_idx], active[b_idx]
                q = (n - 2) * d[(i, j)] - total[i] - total[j]
                if best_q is None or q < best_q:
                    best_q = q
                    best = (i, j)
        i, j = best
        dij = d[(i, j)]
        delta_i = max(0.5 * dij + (total[i] - total[j]) / (2 * (n - 2)), 0.0)
        delta_j = max(dij - delta_i, 0.0)
        counter += 1
        new_node = f"internal_{counter}"
        nodes[new_node] = f"({nodes[i]}:{delta_i:.4f},{nodes[j]}:{delta_j:.4f})"
        for k in active:
            if k in (i, j):
                continue
            d[(new_node, k)] = d[(k, new_node)] = 0.5 * (d[(i, k)] + d[(j, k)] - dij)
        active = [k for k in active if k not in (i, j)] + [new_node]

    i, j = active
    dij = d[(i, j)]
    return f"({nodes[i]}:{dij / 2:.4f},{nodes[j]}:{dij / 2:.4f});"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--distances-json", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = json.loads(Path(args.distances_json).read_text(encoding="utf-8"))
    try:
        newick = neighbor_joining(rows)
        result = {"status": "PASS", "newick": newick}
    except WorkflowError as error:
        result = {"status": "FAIL", "error": str(error)}
    write_json(Path(args.output), result)
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
