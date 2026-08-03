#!/usr/bin/env python3
"""Measure exact globin-target coverage by a benchmark confidence BED."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.truth_regions import (  # noqa: E402
    TruthRegionError,
    audit_regions,
    parse_target,
    read_bed,
    write_audit_tsv,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit 1-based inclusive NAME=chrom:start-end targets against a "
            "0-based half-open benchmark confidence BED."
        )
    )
    parser.add_argument("--bed", required=True)
    parser.add_argument(
        "--target",
        action="append",
        required=True,
        help="Repeatable target specification, e.g. HBA=chr16:170000-178000",
    )
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--summary-json")
    parser.add_argument(
        "--require-full-coverage",
        action="store_true",
        help="Return exit code 3 when any target is not fully covered.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        targets = [parse_target(value) for value in args.target]
        audits = audit_regions(read_bed(args.bed), targets)
        write_audit_tsv(args.output_tsv, audits)
        payload = {
            "schema_version": 1,
            "confidence_bed": str(Path(args.bed)),
            "targets": [item.to_dict() for item in audits],
            "all_fully_covered": all(item.fully_covered for item in audits),
        }
        if args.summary_json:
            summary = Path(args.summary_json)
            summary.parent.mkdir(parents=True, exist_ok=True)
            summary.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
    except (TruthRegionError, OSError, ValueError) as exc:
        print(f"audit_truth_regions: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(payload, sort_keys=True))
    if args.require_full_coverage and not payload["all_fully_covered"]:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
