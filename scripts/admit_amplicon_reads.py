#!/usr/bin/env python3
"""Assign FASTQ molecules to declared/compiled amplicon products."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.assay import AssayProfileError, load_assay_profile  # noqa: E402
from nanoglobin.molecules import (  # noqa: E402
    AdmissionAccumulator,
    AdmissionConfig,
    AdmissionEngine,
    MoleculeAdmissionError,
    OBSERVATION_FIELDS,
    load_compiled_assay,
    read_fastq,
    summary_payload,
    write_product_counts_tsv,
)
import csv  # noqa: E402


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recognise terminal assay primers, orient reads, classify molecule "
            "completeness/artifacts and score complete reads against compiled products."
        )
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fastq", required=True, help="FASTQ or FASTQ.GZ")
    parser.add_argument("--profile", required=True, help="declared assay YAML/JSON")
    parser.add_argument("--compiled-json", required=True)
    parser.add_argument("--molecules-tsv", required=True)
    parser.add_argument("--product-counts-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--end-search-bp", type=int, default=100)
    parser.add_argument("--primer-seed-length", type=int, default=7)
    parser.add_argument("--min-primer-overlap", type=int, default=10)
    parser.add_argument("--max-primer-mismatches", type=int, default=4)
    parser.add_argument("--max-primer-error-rate", type=float, default=0.25)
    parser.add_argument("--length-tolerance-fraction", type=float, default=0.10)
    parser.add_argument("--max-sequence-edit-rate", type=float, default=0.25)
    parser.add_argument("--sequence-assignment-margin", type=float, default=0.02)
    parser.add_argument("--assumed-sequence-error-rate", type=float, default=0.10)
    parser.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="optional deterministic read limit for exploratory runs (0 = all)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = AdmissionConfig(
        end_search_bp=args.end_search_bp,
        primer_seed_length=args.primer_seed_length,
        min_primer_overlap=args.min_primer_overlap,
        max_primer_mismatches=args.max_primer_mismatches,
        max_primer_error_rate=args.max_primer_error_rate,
        length_tolerance_fraction=args.length_tolerance_fraction,
        max_sequence_edit_rate=args.max_sequence_edit_rate,
        sequence_assignment_margin=args.sequence_assignment_margin,
        assumed_sequence_error_rate=args.assumed_sequence_error_rate,
    )
    try:
        config.validate()
        profile = load_assay_profile(args.profile)
        compiled = load_compiled_assay(args.compiled_json, profile)
        engine = AdmissionEngine(profile=profile, compiled=compiled, config=config)
        accumulator = AdmissionAccumulator()
        molecule_path = Path(args.molecules_tsv)
        molecule_path.parent.mkdir(parents=True, exist_ok=True)
        with molecule_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=OBSERVATION_FIELDS, delimiter="\t"
            )
            writer.writeheader()
            for index, record in enumerate(read_fastq(args.fastq), start=1):
                if args.max_reads and index > args.max_reads:
                    break
                observation = engine.admit(record, sample=args.sample)
                writer.writerow(observation.to_dict())
                accumulator.add(observation)
    except (AssayProfileError, MoleculeAdmissionError, OSError, ValueError) as exc:
        print(f"admit_amplicon_reads: {exc}", file=sys.stderr)
        return 2

    write_product_counts_tsv(args.product_counts_tsv, accumulator)
    summary = summary_payload(
        sample=args.sample, profile=profile, observations=accumulator
    )
    summary["admission_config"] = {
        "end_search_bp": config.end_search_bp,
        "primer_seed_length": config.primer_seed_length,
        "min_primer_overlap": config.min_primer_overlap,
        "max_primer_mismatches": config.max_primer_mismatches,
        "max_primer_error_rate": config.max_primer_error_rate,
        "length_tolerance_fraction": config.length_tolerance_fraction,
        "max_sequence_edit_rate": config.max_sequence_edit_rate,
        "sequence_assignment_margin": config.sequence_assignment_margin,
        "assumed_sequence_error_rate": config.assumed_sequence_error_rate,
        "max_reads": args.max_reads,
    }
    summary_path = Path(args.summary_json)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
