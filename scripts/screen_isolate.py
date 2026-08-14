#!/usr/bin/env python3
"""Deterministic end-to-end entrypoint for isolate screening."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from collect_provenance import collect
from classify_with_mlst import classify, route_for_organism
from common import WorkflowError, effective_database_root, write_json
from inspect_input import inspect
from parse_mashpit_results import interpret, load_candidates, locate_candidate_file
from run_assembly_workflow import run_workflow
from run_mashpit import run_mashpit
from validate_assembly import assess


def sample_name(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fasta", ".fna", ".fastq", ".fq", ".fa"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    for mate in ("_R1", "_1", ".R1", ".1", "-R1", "-1"):
        if name.endswith(mate):
            name = name[: -len(mate)]
            break
    return name


def summary_text(result: dict[str, Any]) -> str:
    lines = [f"Sample: {result['sample']}", f"Status: {result['status']}"]
    if result.get("input_type"):
        lines.append(f"Input: {result['input_type']}")
    if result.get("organism"):
        lines.append(f"Organism: {result['organism']}")
    if result.get("assembly_qc"):
        lines.append(f"Assembly QC: {result['assembly_qc']['status']}")
    if result.get("database"):
        lines.append(f"Mashpit database: {result['database']['name']} {result['database']['version']}")
    match = result.get("mashpit_result") or {}
    best = match.get("best_candidate")
    if best:
        lines.extend([
            f"Best candidate cluster: {best['cluster']} (score {best['score']:.3f})",
            f"Screening result: {match['screening_result']}",
            "Interpretation: This isolate is most similar to representatives associated with the candidate cluster.",
            "This screening result does not by itself establish outbreak relatedness.",
            "Recommended next step: Run a validated high-resolution SNP comparison against representative isolates.",
        ])
    if result.get("stop_reason"):
        lines.append(f"Analysis stopped: {result['stop_reason']}")
    warnings = result.get("warnings", [])
    if warnings:
        lines.append("Warnings: " + "; ".join(warnings))
    return "\n".join(lines)


def screen(
    inputs: list[Path], output_dir: Path, database_root: Path,
    organism: str | None = None,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=False)
    commands: list[list[str]] = []
    assembly: Path | None = None
    database_metadata = None
    result: dict[str, Any] = {
        "schema_version": "1.0.0",
        "sample": sample_name(inputs[0]),
        "status": "FAILED",
        "warnings": [],
        "requested_organism": organism,
    }
    try:
        detected = inspect(inputs)
        result["input_type"] = detected["input_type"]
        write_json(output_dir / "input.json", detected)
        if detected["input_type"] == "assembly":
            assembly = Path(detected["assembly_path"])
        else:
            assembly_run = run_workflow(
                Path(detected["r1_path"]), Path(detected["r2_path"]), output_dir / "assembly"
            )
            commands.extend(assembly_run["commands"])
            result["read_qc"] = assembly_run["read_qc"]
            result["warnings"].extend(
                ["Read QC is WARN; downstream interpretation may be affected."]
                if assembly_run["read_qc"]["status"] == "WARN" else []
            )
            assembly = Path(assembly_run["assembly_path"])
        basic_qc = assess(assembly)
        write_json(output_dir / "assembly_qc.initial.json", basic_qc)
        if basic_qc["status"] == "FAIL":
            result["assembly_qc"] = basic_qc
            raise WorkflowError("Assembly QC failed: " + "; ".join(basic_qc["failures"]))
        if result.get("read_qc"):
            total_bases = result["read_qc"].get("total_bases")
            if total_bases:
                coverage = total_bases / basic_qc["metrics"]["total_length"]
                result["read_qc"]["estimated_coverage"] = coverage
                if coverage < 20.0:
                    result["read_qc"]["status"] = "FAIL"
                    raise WorkflowError("Estimated read coverage is below the fixed 20× minimum.")
        if organism:
            routing = route_for_organism(organism)
        else:
            routing = classify(assembly, output_dir / "mlst")
            if routing.get("command"):
                commands.append(routing["command"])
        write_json(output_dir / "routing.json", routing)
        result["routing"] = routing
        if routing["status"] == "UNSUPPORTED":
            raise WorkflowError(
                "The organism does not match one of the currently supported Mashpit databases."
            )
        if routing["status"] != "SUPPORTED":
            raise WorkflowError(
                "Local mlst could not confidently select a supported organism database; "
                "specify --organism only when the organism is known independently."
            )
        result["organism"] = routing["organism"]
        result["warnings"].extend(routing.get("warnings", []))
        routed_qc = assess(assembly, routing["organism_key"])
        write_json(output_dir / "assembly_qc.json", routed_qc)
        result["assembly_qc"] = routed_qc
        if routed_qc["status"] == "FAIL":
            raise WorkflowError("Assembly QC failed: " + "; ".join(routed_qc["failures"]))
        if routed_qc["status"] == "WARN":
            result["warnings"].extend(routed_qc["warnings"])
        mashpit_run = run_mashpit(
            assembly,
            database_root / routing["database_name"],
            routing["database_name"],
            output_dir / "mashpit",
        )
        commands.append(mashpit_run["command"])
        database_metadata = mashpit_run["database"]
        result["database"] = database_metadata
        candidate_file = locate_candidate_file(Path(mashpit_run["output_directory"]))
        parsed = interpret(load_candidates(candidate_file))
        parsed["source_file"] = str(candidate_file)
        write_json(output_dir / "mashpit_result.json", parsed)
        result["mashpit_result"] = parsed
        result["warnings"].extend(parsed.get("warnings", []))
        result["status"] = "COMPLETED_WITH_WARNINGS" if result["warnings"] else "COMPLETED"
    except (WorkflowError, OSError, ValueError) as error:
        result["status"] = "STOPPED"
        result["stop_reason"] = str(error)
    finally:
        result["user_summary"] = summary_text(result)
        write_json(output_dir / "result.json", result)
        provenance = collect(inputs, assembly, database_metadata, commands, result)
        write_json(output_dir / "provenance.json", provenance)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument("--database-root")
    parser.add_argument(
        "--organism",
        choices=("salmonella", "ecoli", "listeria", "campylobacter", "cronobacter"),
        help="Known organism/database key. If omitted, local mlst auto-detects a PubMLST scheme.",
    )
    args = parser.parse_args()
    output_dir = Path(args.output).expanduser().resolve()
    if output_dir.exists():
        print(f"Refusing to overwrite existing output directory: {output_dir}")
        return 2
    database_root = effective_database_root(args.database_root)
    result = screen(
        [Path(value).expanduser().resolve() for value in args.inputs],
        output_dir,
        database_root,
        args.organism,
    )
    print(result["user_summary"])
    return 0 if result["status"].startswith("COMPLETED") else 2


if __name__ == "__main__":
    raise SystemExit(main())
