#!/usr/bin/env python3
"""Validate paired FASTQ records without loading reads into memory."""

from __future__ import annotations

import argparse
import gzip
from contextlib import ExitStack
from pathlib import Path
from typing import TextIO

from common import WorkflowError, write_json


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="strict")
    return path.open(encoding="ascii", errors="strict")


def normalized_read_name(header: str) -> str:
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def validate_pair(r1: Path, r2: Path) -> dict:
    for path in (r1, r2):
        if not path.is_file():
            raise WorkflowError(f"FASTQ does not exist: {path}")
    read_pairs = 0
    total_bases = 0
    errors: list[str] = []
    try:
        with ExitStack() as stack:
            handles = [stack.enter_context(open_text(path)) for path in (r1, r2)]
            while True:
                records = [[handle.readline() for _ in range(4)] for handle in handles]
                if all(not record[0] for record in records):
                    break
                read_pairs += 1
                for mate, record in enumerate(records, 1):
                    if any(line == "" for line in record):
                        errors.append(f"truncated record at pair {read_pairs}, mate {mate}")
                        break
                    header, sequence, plus, quality = [line.rstrip("\r\n") for line in record]
                    if not header.startswith("@"):
                        errors.append(f"invalid header at pair {read_pairs}, mate {mate}")
                    if not plus.startswith("+"):
                        errors.append(f"invalid separator at pair {read_pairs}, mate {mate}")
                    if len(sequence) != len(quality):
                        errors.append(f"sequence/quality length mismatch at pair {read_pairs}, mate {mate}")
                    if not sequence:
                        errors.append(f"empty sequence at pair {read_pairs}, mate {mate}")
                    total_bases += len(sequence)
                if all(record[0] for record in records):
                    names = [normalized_read_name(record[0].rstrip()) for record in records]
                    if names[0] != names[1]:
                        errors.append(f"mate names differ at pair {read_pairs}")
                if errors:
                    break
    except (OSError, UnicodeError, EOFError) as error:
        errors.append(f"cannot decode FASTQ: {error}")
    if read_pairs == 0:
        errors.append("FASTQ pair contains no reads")
    return {
        "status": "FAIL" if errors else "PASS",
        "read_pairs": read_pairs,
        "total_bases": total_bases,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("r1")
    parser.add_argument("r2")
    parser.add_argument("--output")
    args = parser.parse_args()
    try:
        result = validate_pair(Path(args.r1), Path(args.r2))
    except WorkflowError as error:
        result = {"status": "FAIL", "read_pairs": 0, "total_bases": 0, "errors": [str(error)]}
    if args.output:
        write_json(Path(args.output), result)
    else:
        import json

        print(json.dumps(result, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
