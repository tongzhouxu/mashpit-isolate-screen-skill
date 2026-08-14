#!/usr/bin/env python3
"""Detect one assembly or one paired-end Illumina FASTQ pair."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from common import WorkflowError, write_json


ASSEMBLY_SUFFIXES = (".fa", ".fasta", ".fna")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
MATE_PATTERN = re.compile(r"(?i)(.*?)(?:[_\.-]R?)([12])(?:[_\.-]?\d+)?(\.f(?:ast)?q(?:\.gz)?)$")


def inspect(paths: list[Path]) -> dict:
    resolved = [path.expanduser().resolve() for path in paths]
    missing = [str(path) for path in resolved if not path.is_file()]
    if missing:
        raise WorkflowError("Input file does not exist: " + ", ".join(missing))
    if len(resolved) == 1 and resolved[0].name.lower().endswith(ASSEMBLY_SUFFIXES):
        return {"input_type": "assembly", "assembly_path": str(resolved[0])}
    if len(resolved) != 2 or not all(path.name.lower().endswith(FASTQ_SUFFIXES) for path in resolved):
        raise WorkflowError("Provide exactly one FASTA/FNA assembly or exactly two paired FASTQ files.")
    matches = [MATE_PATTERN.match(path.name) for path in resolved]
    if not all(matches):
        raise WorkflowError("FASTQ names must identify mates with R1/R2 or _1/_2.")
    parsed = [(match.group(1), int(match.group(2))) for match in matches if match]
    if parsed[0][0].lower() != parsed[1][0].lower() or {parsed[0][1], parsed[1][1]} != {1, 2}:
        raise WorkflowError("FASTQ files do not appear to be a matching R1/R2 pair.")
    by_mate = {mate: resolved[index] for index, (_, mate) in enumerate(parsed)}
    return {
        "input_type": "illumina_paired_fastq",
        "r1_path": str(by_mate[1]),
        "r2_path": str(by_mate[2]),
        "sample": parsed[0][0].rstrip("_.-"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--output")
    args = parser.parse_args()
    try:
        result = inspect([Path(value) for value in args.inputs])
    except WorkflowError as error:
        result = {"status": "FAIL", "error": str(error)}
        if args.output:
            write_json(Path(args.output), result)
        print(str(error))
        return 2
    if args.output:
        write_json(Path(args.output), result)
    else:
        import json

        print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
