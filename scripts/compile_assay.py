#!/usr/bin/env python3
"""Compile a declared amplicon assay against chromosome-haplotype sequences."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Allow execution from a repository checkout without installation.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.assay import (  # noqa: E402
    AssayProfileError,
    compile_catalogue,
    coverage_rows,
    indistinguishability_classes,
    load_assay_profile,
    read_fasta,
    write_compiled_json,
    write_coverage_tsv,
    write_indistinguishability_tsv,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run in-silico PCR for every declared assay product and haplotype, "
            "then report assay coverage and diploid genotype ambiguity classes."
        )
    )
    parser.add_argument("--profile", required=True, help="YAML/JSON assay profile")
    parser.add_argument("--haplotypes", required=True, help="chromosome-haplotype FASTA")
    parser.add_argument("--compiled-json", required=True, help="compiled product JSON output")
    parser.add_argument("--coverage-tsv", required=True, help="haplotype/product coverage TSV")
    parser.add_argument(
        "--indistinguishability-tsv",
        required=True,
        help="diploid genotype indistinguishability classes",
    )
    parser.add_argument(
        "--summary-json",
        help="optional compact summary for workflow logs and regression checks",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        profile = load_assay_profile(args.profile)
        haplotypes = read_fasta(args.haplotypes)
        compiled = compile_catalogue(profile, haplotypes)
        coverage = list(coverage_rows(profile, haplotypes, compiled))
        ambiguity = indistinguishability_classes(haplotypes, compiled)
    except (AssayProfileError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"compile_assay: {exc}", file=sys.stderr)
        return 2

    write_compiled_json(args.compiled_json, profile, haplotypes, compiled)
    write_coverage_tsv(args.coverage_tsv, coverage)
    write_indistinguishability_tsv(args.indistinguishability_tsv, ambiguity)

    required_missing = sum(
        1 for row in coverage if row["required"] and not row["observable"]
    )
    summary = {
        "assay_id": profile.assay_id,
        "assay_version": profile.version,
        "haplotypes": len(haplotypes),
        "product_definitions": len(profile.products),
        "compiled_products": len(compiled),
        "required_product_gaps": required_missing,
        "indistinguishability_classes": len(ambiguity),
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
