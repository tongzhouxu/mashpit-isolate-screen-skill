#!/usr/bin/env python3
"""Render the ska2 Neighbor-Joining tree to PNG, with the query genome highlighted.

Uses phytreeviz + matplotlib, both already pulled in transitively by the
pinned mashpit commit's own dependencies (phytreeviz ~= 0.2.0 in mashpit's
pyproject.toml, which itself requires matplotlib) - no new container pin
needed. Rendering is best-effort: a failure here degrades to a warning, it
never fails the underlying SNP resolution or Mash screen.
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

from common import write_json


QUERY_SAMPLE_ID = "QUERY"
HIGHLIGHT_COLOR = "#d62728"


def render(newick: str, output_path: Path, taxa_count: int) -> None:
    from phytreeviz import TreeViz

    # Scale down per-leaf spacing/label size as the tree grows so large trees
    # (up to max_total_genomes + 1 taxa) stay readable rather than enormous.
    height = 0.45 if taxa_count <= 30 else max(0.12, 12.0 / taxa_count)
    label_size = 9 if taxa_count <= 30 else max(4, int(220 / taxa_count))

    with tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as handle:
        handle.write(newick)
        nwk_path = Path(handle.name)
    try:
        tree_viz = TreeViz(str(nwk_path), height=height, width=8, leaf_label_size=label_size)
        tree_viz.highlight(QUERY_SAMPLE_ID, HIGHLIGHT_COLOR, area="branch-label")
        tree_viz.set_node_label_props(QUERY_SAMPLE_ID, color=HIGHLIGHT_COLOR, fontweight="bold")
        tree_viz.show_scale_bar()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        tree_viz.savefig(str(output_path))
    finally:
        nwk_path.unlink(missing_ok=True)


def render_from_interpretation(interpretation: dict, output_path: Path) -> dict:
    newick = interpretation.get("newick_tree")
    if not newick:
        return {"status": "SKIPPED", "reason": "No newick_tree in the SNP interpretation."}
    taxa_count = len(interpretation.get("ranked", [])) + 1
    try:
        render(newick, output_path, taxa_count)
        return {"status": "PASS", "path": str(output_path)}
    except Exception as error:  # rendering is best-effort; never raise into the caller
        return {"status": "FAIL", "error": str(error)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--interpretation-json", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    interpretation = json.loads(Path(args.interpretation_json).read_text(encoding="utf-8"))
    result = render_from_interpretation(interpretation, Path(args.output))
    write_json(Path(args.output).with_suffix(".render.json"), result)
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
