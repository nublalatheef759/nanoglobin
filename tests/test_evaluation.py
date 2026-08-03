from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

import scripts.check_evaluation as evaluation
from scripts.check_evaluation import (
    EvaluationError,
    summarize,
    validate_cython_proof,
    validate_happy_summary,
    validate_public_runs,
    validate_small_variants,
    validate_structural_comparators,
)


def test_real_small_variant_snapshots_are_arithmetically_consistent() -> None:
    rows = validate_small_variants()
    by_id = {row["evaluation_id"]: row for row in rows}
    core = by_id["pooled_core_snv"]
    assert core["tp"] == 140
    assert core["fn"] == 5
    assert core["fp"] == 1
    assert math.isclose(core["recall"], 140 / 145)
    assert math.isclose(core["precision"], 140 / 141)
    assert set(core["cohort"].split(";")) == {
        "HG002", "HG02071", "HG02083", "HG02514", "HG02074", "HG02622"
    }
    giab = by_id["giab_hg002_globin_small_variants"]
    assert giab["cohort"] == "HG002"
    assert (giab["tp"], giab["fn"], giab["fp"]) == (23, 0, 0)


def test_happy_source_snapshot_matches_the_committed_counts() -> None:
    rows = validate_small_variants()
    summary = validate_happy_summary(rows)
    assert "Pooled across HG002 + 5 HPRC genomes." in summary
    assert "140 TP, 5 FN, 1 FP" in summary
    assert "31 TP, 18 FN, 3 FP" in summary


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_happy_snapshot_rejects_arithmetically_valid_structured_drift(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = evaluation.read_tsv(evaluation.SMALL_VARIANTS)
    changed = [dict(row) for row in rows]
    core = next(
        row for row in changed if row["evaluation_id"] == "pooled_core_snv"
    )
    core.update(
        tp="141",
        fn="4",
        fp="1",
        recall=str(141 / 145),
        precision=str(141 / 142),
    )
    changed_path = tmp_path / "small_variant_results.tsv"
    _write_rows(changed_path, changed)
    monkeypatch.setattr(evaluation, "SMALL_VARIANTS", changed_path)

    validated = evaluation.validate_small_variants()
    with pytest.raises(
        EvaluationError,
        match=r"pooled_core_snv tp: snapshot=140, structured=141",
    ):
        evaluation.validate_happy_summary(validated)


def test_happy_snapshot_rejects_source_text_drift(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    changed_summary = tmp_path / "happy_summary_final.txt"
    text = evaluation.HAPPY_SUMMARY.read_text(encoding="utf-8")
    changed_summary.write_text(
        text.replace(
            "SNV:   96.6% recall, 99.3% precision (140 TP, 5 FN, 1 FP)",
            "SNV:   96.6% recall, 99.3% precision (139 TP, 6 FN, 1 FP)",
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(evaluation, "HAPPY_SUMMARY", changed_summary)

    with pytest.raises(
        EvaluationError,
        match=r"pooled_core_snv tp: snapshot=139, structured=140",
    ):
        evaluation.validate_happy_summary(evaluation.validate_small_variants())


def test_structural_rows_remain_comparator_only() -> None:
    rows = validate_structural_comparators()
    assert {row["truth_role"] for row in rows} == {"comparator_only"}
    assert {row["sample_id"] for row in rows} >= {
        "NA21106",
        "HG00642",
        "HG03136",
        "HG00735",
        "HG03862",
    }


def test_public_amplicon_metadata_is_real_and_unique() -> None:
    rows = validate_public_runs()
    assert len(rows) == 64
    assert {row["instrument_model"] for row in rows} == {"GridION"}
    assert len({row["run_accession"] for row in rows}) == 64


def test_cython_proof_is_internally_consistent() -> None:
    proof = validate_cython_proof()
    assert proof["workload"]["candidate_haplotypes"] == 24
    assert proof["workload"]["observable_genotype_classes"] == 300
    assert proof["workload"]["molecule_rows"] == 12000
    assert proof["median_speedup"] > 100


def test_repository_evaluation_summary() -> None:
    summary = summarize()
    assert summary["status"] == "ok"
    assert summary["small_variant_evaluations"] == 6
    assert len(summary["small_variant_samples"]) == 6
    assert summary["structural_comparator_rows"] == 6
    assert summary["public_amplicon_runs"] == 64
