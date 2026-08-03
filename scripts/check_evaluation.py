#!/usr/bin/env python3
"""Validate and summarize committed, non-synthetic evaluation snapshots.

This command does not rerun hap.py or download the original WGS data. It verifies
that the committed aggregate counts are internally consistent, that the structured
small-variant table agrees with the retained hap.py source snapshot, that
comparator-only HBA records are not promoted to truth, and that the retained Cython
benchmark proof matches its raw timing fields.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
SMALL_VARIANTS = ROOT / "evaluation" / "small_variant_results.tsv"
HAPPY_SUMMARY = ROOT / "evaluation" / "happy_summary_final.txt"
STRUCTURAL = ROOT / "evaluation" / "hba_structural_comparators.tsv"
PUBLIC_RUNS = ROOT / "evaluation" / "public_amplicon_runs.tsv"
CYTHON_PROOF = ROOT / "benchmarks" / "cython_genotype_proof.json"

_HAPPY_SECTION_PREFIXES = {
    "CORE (HBA/HBB):": "pooled_core",
    (
        "EXTENDED beta-cluster, GENE BODIES "
        "(HBB,HBD,HBG1/2,HBE1 - clinically relevant):"
    ): "pooled_beta_gene_body",
    (
        "EXTENDED beta-cluster, FULL interval "
        "(incl. intergenic + LCR):"
    ): "pooled_beta_full_interval",
}
_HAPPY_RESULT = re.compile(
    r"^(?P<variant_type>SNV|INDEL):\s+"
    r"(?P<recall>[0-9]+(?:\.[0-9]+)?)% recall,\s+"
    r"(?P<precision>[0-9]+(?:\.[0-9]+)?)% precision"
    r"(?:\s+\((?P<tp>[0-9]+) TP,\s*(?P<fn>[0-9]+) FN,\s*"
    r"(?P<fp>[0-9]+) FP\))?$"
)
_EXPECTED_HAPPY_EVALUATIONS = {
    "pooled_core_snv",
    "pooled_core_indel",
    "pooled_beta_gene_body_snv",
    "pooled_beta_gene_body_indel",
    "pooled_beta_full_interval_snv",
}


class EvaluationError(ValueError):
    """Raised when a committed evaluation snapshot is inconsistent."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        try:
            label = str(path.relative_to(ROOT))
        except ValueError:
            label = str(path)
        raise EvaluationError(f"{label} contains no data rows")
    return rows


def require_unique(rows: Iterable[dict[str, str]], field: str, label: str) -> None:
    values = [row.get(field, "").strip() for row in rows]
    if any(not value for value in values):
        raise EvaluationError(f"{label} contains an empty {field}")
    if len(values) != len(set(values)):
        raise EvaluationError(f"{label} contains duplicate {field} values")


def close(left: float, right: float, tolerance: float = 1e-12) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def validate_small_variants() -> list[dict[str, Any]]:
    rows = read_tsv(SMALL_VARIANTS)
    require_unique(rows, "evaluation_id", "small-variant results")
    output: list[dict[str, Any]] = []
    for line_number, row in enumerate(rows, start=2):
        try:
            tp = int(row["tp"])
            fn = int(row["fn"])
            fp = int(row["fp"])
            recall = float(row["recall"])
            precision = float(row["precision"])
        except (KeyError, ValueError) as exc:
            raise EvaluationError(
                f"invalid numeric field in small-variant results:{line_number}"
            ) from exc
        if min(tp, fn, fp) < 0:
            raise EvaluationError(f"negative count at small-variant row {line_number}")
        source_script = row.get("source_script", "").strip()
        if not source_script or not (ROOT / source_script).is_file():
            raise EvaluationError(
                f"missing source script for {row.get('evaluation_id', line_number)!r}: "
                f"{source_script!r}"
            )
        if not row.get("truth_source", "").strip():
            raise EvaluationError(
                f"missing truth source at small-variant row {line_number}"
            )
        if not row.get("comparison", "").strip():
            raise EvaluationError(
                f"missing comparison method at small-variant row {line_number}"
            )
        expected_recall = tp / (tp + fn) if tp + fn else float("nan")
        expected_precision = tp / (tp + fp) if tp + fp else float("nan")
        if not close(recall, expected_recall):
            raise EvaluationError(
                f"recall mismatch for {row['evaluation_id']}: {recall} != {expected_recall}"
            )
        if not close(precision, expected_precision):
            raise EvaluationError(
                f"precision mismatch for {row['evaluation_id']}: "
                f"{precision} != {expected_precision}"
            )
        output.append(
            {
                **row,
                "tp": tp,
                "fn": fn,
                "fp": fp,
                "recall": recall,
                "precision": precision,
            }
        )
    return output


def parse_happy_summary(text: str) -> dict[str, dict[str, int | float]]:
    """Parse count-bearing hap.py records from the retained human-readable snapshot."""

    if "Pooled across HG002 + 5 HPRC genomes." not in text:
        raise EvaluationError(
            "hap.py source snapshot is missing its pooled-cohort declaration"
        )

    current_prefix: str | None = None
    records: dict[str, dict[str, int | float]] = {}
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if line in _HAPPY_SECTION_PREFIXES:
            current_prefix = _HAPPY_SECTION_PREFIXES[line]
            continue

        match = _HAPPY_RESULT.fullmatch(line)
        if match is None:
            continue
        if current_prefix is None:
            raise EvaluationError(
                f"hap.py result outside a recognised section at line {line_number}"
            )
        if match.group("tp") is None:
            # The retained full-interval INDEL line contains rounded percentages
            # only, so it cannot be the source of a structured TP/FN/FP record.
            continue

        evaluation_id = (
            f"{current_prefix}_{match.group('variant_type').lower()}"
        )
        if evaluation_id in records:
            raise EvaluationError(
                f"duplicate hap.py source record for {evaluation_id}"
            )
        records[evaluation_id] = {
            "tp": int(match.group("tp")),
            "fn": int(match.group("fn")),
            "fp": int(match.group("fp")),
            "recall_percent": float(match.group("recall")),
            "precision_percent": float(match.group("precision")),
        }

    found = set(records)
    if found != _EXPECTED_HAPPY_EVALUATIONS:
        missing = sorted(_EXPECTED_HAPPY_EVALUATIONS - found)
        extra = sorted(found - _EXPECTED_HAPPY_EVALUATIONS)
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unexpected " + ", ".join(extra))
        raise EvaluationError(
            "hap.py source snapshot records are incomplete: " + "; ".join(details)
        )
    return records


def _cross_check_happy_rows(
    source_records: Mapping[str, Mapping[str, int | float]],
    small_variants: Sequence[Mapping[str, Any]],
) -> None:
    structured = {
        str(row["evaluation_id"]): row
        for row in small_variants
        if str(row.get("source_script", "")).strip() == "run_happy.sh"
    }
    source_ids = set(source_records)
    structured_ids = set(structured)
    if source_ids != structured_ids:
        missing = sorted(source_ids - structured_ids)
        extra = sorted(structured_ids - source_ids)
        details = []
        if missing:
            details.append("missing structured rows: " + ", ".join(missing))
        if extra:
            details.append("rows without source counts: " + ", ".join(extra))
        raise EvaluationError(
            "hap.py snapshot/structured-row coverage mismatch: "
            + "; ".join(details)
        )

    for evaluation_id in sorted(source_records):
        source = source_records[evaluation_id]
        row = structured[evaluation_id]
        for field in ("tp", "fn", "fp"):
            source_value = int(source[field])
            row_value = int(row[field])
            if source_value != row_value:
                raise EvaluationError(
                    f"hap.py source snapshot mismatch for {evaluation_id} {field}: "
                    f"snapshot={source_value}, structured={row_value}"
                )

        expected_recall_percent = round(100.0 * float(row["recall"]), 1)
        expected_precision_percent = round(100.0 * float(row["precision"]), 1)
        if not close(
            float(source["recall_percent"]), expected_recall_percent
        ):
            raise EvaluationError(
                f"hap.py source snapshot mismatch for {evaluation_id} recall: "
                f"snapshot={source['recall_percent']}%, "
                f"structured={expected_recall_percent}%"
            )
        if not close(
            float(source["precision_percent"]), expected_precision_percent
        ):
            raise EvaluationError(
                f"hap.py source snapshot mismatch for {evaluation_id} precision: "
                f"snapshot={source['precision_percent']}%, "
                f"structured={expected_precision_percent}%"
            )


def validate_happy_summary(
    small_variants: Sequence[Mapping[str, Any]] | None = None,
) -> str:
    text = HAPPY_SUMMARY.read_text(encoding="utf-8")
    source_records = parse_happy_summary(text)
    _cross_check_happy_rows(
        source_records,
        list(small_variants) if small_variants is not None else validate_small_variants(),
    )
    return text


def validate_structural_comparators() -> list[dict[str, str]]:
    rows = read_tsv(STRUCTURAL)
    require_unique(rows, "sample_id", "HBA structural comparator results")
    for line_number, row in enumerate(rows, start=2):
        if row.get("truth_role") != "comparator_only":
            raise EvaluationError(
                f"structural row {line_number} must remain comparator_only"
            )
        event = row.get("observed_cigar_event_bp", "").strip()
        fraction = row.get("deletion_read_fraction", "").strip()
        if event and int(event) <= 0:
            raise EvaluationError(f"invalid event size at structural row {line_number}")
        if fraction:
            value = float(fraction)
            if not 0 <= value <= 1:
                raise EvaluationError(
                    f"invalid deletion fraction at structural row {line_number}"
                )
    return rows


def validate_public_runs() -> list[dict[str, str]]:
    rows = read_tsv(PUBLIC_RUNS)
    require_unique(rows, "run_accession", "public amplicon run metadata")
    for line_number, row in enumerate(rows, start=2):
        if not row.get("instrument_model", "").strip():
            raise EvaluationError(
                f"missing instrument model in public-run row {line_number}"
            )
    return rows


def validate_cython_proof() -> dict[str, Any]:
    proof = json.loads(CYTHON_PROOF.read_text(encoding="utf-8"))
    if not str(proof.get("numerical_parity", "")).startswith("passed"):
        raise EvaluationError("Cython proof does not record numerical parity")
    python_median = float(proof["historical_python_seconds"]["median"])
    cython_median = float(proof["cython_seconds"]["median"])
    speedup = float(proof["median_speedup"])
    expected = python_median / cython_median
    if not math.isclose(speedup, expected, rel_tol=1e-12, abs_tol=1e-12):
        raise EvaluationError(
            f"Cython speedup mismatch: {speedup} != {expected}"
        )
    if speedup < 2:
        raise EvaluationError("Cython replacement proof is below the 2x gate")
    return proof


def summarize() -> dict[str, Any]:
    small_variants = validate_small_variants()
    validate_happy_summary(small_variants)
    structural = validate_structural_comparators()
    public_runs = validate_public_runs()
    proof = validate_cython_proof()
    return {
        "schema_version": 1,
        "small_variant_evaluations": len(small_variants),
        "small_variant_samples": sorted(
            {
                sample
                for row in small_variants
                for sample in row["cohort"].split(";")
                if sample
            }
        ),
        "structural_comparator_rows": len(structural),
        "public_amplicon_runs": len(public_runs),
        "public_instruments": sorted(
            {row["instrument_model"] for row in public_runs}
        ),
        "cython_workload": proof["workload"],
        "cython_median_speedup": proof["median_speedup"],
        "status": "ok",
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate committed real-evaluation snapshots and provenance roles."
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit only machine-readable JSON",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        summary = summarize()
    except (EvaluationError, OSError, KeyError, ValueError, json.JSONDecodeError) as exc:
        print(f"check_evaluation: {exc}")
        return 2
    if args.json:
        print(json.dumps(summary, sort_keys=True))
    else:
        print("evaluation snapshots: OK")
        print(f"  small-variant evaluations: {summary['small_variant_evaluations']}")
        print(f"  samples represented: {len(summary['small_variant_samples'])}")
        print(f"  structural comparator rows: {summary['structural_comparator_rows']}")
        print(f"  public amplicon runs: {summary['public_amplicon_runs']}")
        print(f"  Cython median speedup gate: {summary['cython_median_speedup']:.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
