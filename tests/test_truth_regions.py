from __future__ import annotations

import json
from pathlib import Path

import pytest

from nanoglobin.truth_regions import (
    TruthRegionError,
    audit_regions,
    parse_target,
    read_bed,
)
from scripts.audit_truth_regions import main as audit_main


def test_merges_confidence_intervals_and_measures_partial_coverage(
    tmp_path: Path,
) -> None:
    bed = tmp_path / "confidence.bed"
    bed.write_text(
        "chr16\t99\t150\n"
        "chr16\t140\t180\n"
        "chr16\t200\t210\n"
        "chr11\t10\t20\n",
        encoding="utf-8",
    )
    confidence = read_bed(bed)
    audits = audit_regions(
        confidence,
        [
            parse_target("HBA=chr16:100-200"),
            parse_target("HBB=chr11:11-20"),
        ],
    )
    hba, hbb = audits
    assert hba.length_bp == 101
    assert hba.covered_bp == 81
    assert hba.uncovered_bp == 20
    assert not hba.fully_covered
    assert hbb.covered_bp == 10
    assert hbb.fully_covered


def test_invalid_or_duplicate_targets_are_rejected(tmp_path: Path) -> None:
    with pytest.raises(TruthRegionError, match="expected NAME"):
        parse_target("chr16:1-2")
    target = parse_target("HBA=chr16:1-2")
    with pytest.raises(TruthRegionError, match="duplicate target"):
        audit_regions({"chr16": [(0, 2)]}, [target, target])


def test_cli_can_fail_closed_on_incomplete_truth_region(tmp_path: Path) -> None:
    bed = tmp_path / "confidence.bed"
    output = tmp_path / "audit.tsv"
    summary = tmp_path / "summary.json"
    bed.write_text("chr16\t100\t150\n", encoding="utf-8")
    code = audit_main(
        [
            "--bed",
            str(bed),
            "--target",
            "HBA=chr16:101-200",
            "--output-tsv",
            str(output),
            "--summary-json",
            str(summary),
            "--require-full-coverage",
        ]
    )
    assert code == 3
    assert not json.loads(summary.read_text(encoding="utf-8"))[
        "all_fully_covered"
    ]
    assert output.read_text(encoding="utf-8").startswith("name\tchrom")
