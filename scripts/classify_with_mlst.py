#!/usr/bin/env python3
"""Run local MLST scheme detection and map it to a supported Mashpit database."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from common import CONFIG_DIR, WorkflowError, load_json, require_executable, run_logged, write_json


def database_routes() -> dict[str, dict[str, str]]:
    workflow = load_json(CONFIG_DIR / "workflow.json")
    return {item["organism_key"]: item for item in workflow["databases"]}


def route_for_organism(organism_key: str, source: str = "user") -> dict[str, Any]:
    route = database_routes().get(organism_key)
    if route is None:
        raise WorkflowError(f"Unsupported organism key: {organism_key}")
    return {
        "status": "SUPPORTED",
        "source": source,
        "organism_key": route["organism_key"],
        "organism": route["organism"],
        "database_name": route["name"],
        "warnings": [],
    }


def parse_mlst_csv(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise WorkflowError("Expected exactly one result row from local mlst classification.")
    row = {str(key).upper(): (value or "").strip() for key, value in rows[0].items()}
    required = {"FILE", "SCHEME", "ST", "STATUS", "SCORE", "ALLELES"}
    if not required.issubset(row):
        raise WorkflowError("Local mlst output does not match the pinned --full --csv schema.")
    return row


def interpret_mlst(row: dict[str, str]) -> dict[str, Any]:
    policy = load_json(CONFIG_DIR / "mlst-routing.json")
    scheme = row["SCHEME"]
    native_status = row["STATUS"].upper()
    base = {
        "source": "mlst",
        "mlst": {
            "scheme": scheme,
            "sequence_type": row["ST"],
            "status": native_status,
            "score": row["SCORE"],
            "alleles": row["ALLELES"],
        },
        "warnings": [],
    }
    if native_status in policy["rejected_statuses"] or scheme in {"", "-"}:
        return {**base, "status": "UNCERTAIN"}
    organism_key = policy["schemes"].get(scheme)
    if organism_key is None:
        return {**base, "status": "UNSUPPORTED"}
    if native_status not in policy["accepted_statuses"] + policy["warning_statuses"]:
        return {**base, "status": "UNCERTAIN"}
    result = {**route_for_organism(organism_key, "mlst"), **base}
    result["status"] = "SUPPORTED"
    if native_status in policy["warning_statuses"]:
        result["warnings"] = [
            f"mlst selected scheme '{scheme}' with native status {native_status}; database routing may be less reliable."
        ]
    return result


def classify(assembly: Path, output_dir: Path) -> dict[str, Any]:
    if not assembly.is_file():
        raise WorkflowError(f"Assembly is missing: {assembly}")
    mlst = require_executable("mlst")
    output_dir.mkdir(parents=True, exist_ok=True)
    result_path = output_dir / "mlst.csv"
    command = [mlst, "--full", "--csv", "--quiet", str(assembly)]
    returncode = run_logged(command, output_dir, result_path, output_dir / "mlst.stderr.log")
    if returncode:
        raise WorkflowError(f"Local mlst classification failed; inspect {output_dir / 'mlst.stderr.log'}.")
    result = interpret_mlst(parse_mlst_csv(result_path))
    result["command"] = command
    result["output_file"] = str(result_path)
    write_json(output_dir / "mlst_result.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("assembly")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        result = classify(Path(args.assembly), Path(args.output))
    except (WorkflowError, OSError, ValueError) as error:
        result = {"status": "ERROR", "error": str(error)}
        write_json(Path(args.output) / "mlst_result.json", result)
    print(result["status"])
    return 0 if result["status"] == "SUPPORTED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
