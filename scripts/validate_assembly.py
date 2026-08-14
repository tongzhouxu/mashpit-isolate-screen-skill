#!/usr/bin/env python3
"""Validate FASTA and calculate deterministic assembly QC metrics."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import CONFIG_DIR, WorkflowError, load_json, write_json


VALID_BASES = frozenset("ACGTURYSWKMBDHVNacgturyswkmbdhvn-.")


def fasta_lengths(path: Path) -> tuple[list[int], int, list[str]]:
    lengths: list[int] = []
    current = 0
    ambiguous = 0
    errors: list[str] = []
    saw_header = False
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if saw_header:
                    if current == 0:
                        errors.append(f"empty sequence before line {line_number}")
                    lengths.append(current)
                if not line[1:].strip():
                    errors.append(f"empty FASTA header at line {line_number}")
                saw_header = True
                current = 0
                continue
            if not saw_header:
                errors.append(f"sequence occurs before first header at line {line_number}")
                break
            invalid = set(line) - VALID_BASES
            if invalid:
                errors.append(f"invalid sequence symbols at line {line_number}: {''.join(sorted(invalid))}")
                continue
            current += len(line)
            ambiguous += sum(base.upper() not in {"A", "C", "G", "T"} for base in line)
    if saw_header:
        if current == 0:
            errors.append("final FASTA record has no sequence")
        lengths.append(current)
    if not saw_header:
        errors.append("no FASTA records found")
    return lengths, ambiguous, errors


def n50(lengths: list[int]) -> int:
    total = sum(lengths)
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running * 2 >= total:
            return length
    return 0


def assess(path: Path, organism_key: str | None = None) -> dict:
    if not path.is_file():
        raise WorkflowError(f"Assembly does not exist: {path}")
    try:
        lengths, ambiguous, format_errors = fasta_lengths(path)
    except (OSError, UnicodeError) as error:
        raise WorkflowError(f"Cannot read assembly as text FASTA: {error}") from error
    total = sum(lengths)
    metrics = {
        "total_length": total,
        "contig_count": len(lengths),
        "n50": n50(lengths),
        "largest_contig": max(lengths, default=0),
        "ambiguous_bases": ambiguous,
        "ambiguous_fraction": ambiguous / total if total else 1.0,
    }
    policy = load_json(CONFIG_DIR / "qc-policy.json")
    generic = policy["generic"]
    failures = list(format_errors)
    warnings: list[str] = []
    if total < generic["minimum_total_length_fail"] or total > generic["maximum_total_length_fail"]:
        failures.append("assembly length is outside the broad bacterial-isolate safety envelope")
    if metrics["n50"] < generic["minimum_n50_fail"]:
        failures.append(f"N50 is below {generic['minimum_n50_fail']:,} bp")
    elif metrics["n50"] < generic["minimum_n50_warn"]:
        warnings.append(f"N50 is below {generic['minimum_n50_warn']:,} bp")
    max_contigs = generic["maximum_contigs_fail"]
    if organism_key:
        organism_policy = policy["organisms"].get(organism_key)
        if not organism_policy:
            failures.append(f"no QC policy exists for routed organism '{organism_key}'")
        else:
            max_contigs = organism_policy["max_contigs"]
            if not organism_policy["minimum_length"] <= total <= organism_policy["maximum_length"]:
                failures.append("assembly length is outside the routed organism range")
    if metrics["contig_count"] > max_contigs:
        failures.append(f"contig count exceeds {max_contigs}")
    elif metrics["contig_count"] > generic["maximum_contigs_warn"]:
        warnings.append(f"contig count exceeds {generic['maximum_contigs_warn']}")
    if metrics["ambiguous_fraction"] > generic["maximum_ambiguous_fraction_fail"]:
        failures.append("ambiguous-base fraction exceeds 5%")
    elif metrics["ambiguous_fraction"] > generic["maximum_ambiguous_fraction_warn"]:
        warnings.append("ambiguous-base fraction exceeds 1%")
    status = "FAIL" if failures else "WARN" if warnings else "PASS"
    return {
        "status": status,
        "policy_version": policy["policy_version"],
        "metrics": metrics,
        "failures": failures,
        "warnings": warnings,
        "organism_key": organism_key,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("assembly")
    parser.add_argument("--organism-key")
    parser.add_argument("--output")
    args = parser.parse_args()
    try:
        result = assess(Path(args.assembly), args.organism_key)
    except WorkflowError as error:
        result = {"status": "FAIL", "failures": [str(error)], "warnings": [], "metrics": {}}
    if args.output:
        write_json(Path(args.output), result)
    else:
        import json

        print(json.dumps(result, indent=2))
    return 0 if result["status"] != "FAIL" else 2


if __name__ == "__main__":
    raise SystemExit(main())
