#!/usr/bin/env python3
"""Split a wide report-derived cohort CSV into auditable long-form tables."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.cohort import (  # noqa: E402
    CohortContractError,
    load_cohort_contract,
    normalise_rows,
    read_source_csv,
    write_records,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Normalise a selectively reported clinical spreadsheet without "
            "mistaking missing text for a reference genotype."
        )
    )
    parser.add_argument("--input", required=True, help="source CSV export")
    parser.add_argument("--contract", required=True, help="YAML/JSON column contract")
    parser.add_argument("--episodes", required=True, help="normalised episode CSV")
    parser.add_argument("--findings", required=True, help="normalised finding CSV")
    parser.add_argument("--phenotypes", required=True, help="long phenotype CSV")
    parser.add_argument(
        "--source-label",
        help="stable, non-identifying source label; defaults to the input filename",
    )
    parser.add_argument("--summary-json", help="optional audit summary")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    source_label = args.source_label or Path(args.input).name
    try:
        contract = load_cohort_contract(args.contract)
        rows = read_source_csv(args.input, contract)
        episodes, findings, phenotypes = normalise_rows(
            rows, contract, source_label=source_label
        )
    except (CohortContractError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"normalise_cohort_reports: {exc}", file=sys.stderr)
        return 2

    write_records(args.episodes, episodes)
    write_records(args.findings, findings)
    write_records(args.phenotypes, phenotypes)
    summary = {
        "contract_id": contract.contract_id,
        "source_label": source_label,
        "source_rows": len(rows),
        "episodes": len(episodes),
        "findings": len(findings),
        "phenotypes": len(phenotypes),
        "episodes_without_person_id": sum(
            1 for episode in episodes if not episode.get("person_id")
        ),
        "not_tested_findings": sum(
            1 for finding in findings if finding["reporting_status"] == "not_tested"
        ),
        "not_reported_findings": sum(
            1 for finding in findings if finding["reporting_status"] == "not_reported"
        ),
    }
    if args.summary_json:
        Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary_json).write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
