from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from nanoglobin.genotype import (
    GenotypeConfig,
    load_compiled_space,
    make_call,
    read_molecule_evidence,
    score_genotypes,
)
from scripts.genotype_amplicons import main as genotype_main


def write_compiled(path: Path, haplotype_products: dict[str, list[tuple[str, str]]]) -> None:
    products = []
    for haplotype_id, keys in haplotype_products.items():
        for product_id, sequence_hash in keys:
            products.append(
                {
                    "assay_id": "test",
                    "assay_version": "1",
                    "product_id": product_id,
                    "pool": "alpha",
                    "roles": ["whole_gene"],
                    "required": True,
                    "haplotype_id": haplotype_id,
                    "start0": 0,
                    "end0": 100,
                    "length_bp": 100,
                    "forward_mismatches": 0,
                    "reverse_mismatches": 0,
                    "sequence": "A" * 100,
                    "sequence_sha256": sequence_hash,
                }
            )
    payload = {
        "schema_version": 1,
        "assay": {"assay_id": "test", "version": "1", "platform": "ont"},
        "haplotypes": [
            {
                "haplotype_id": haplotype_id,
                "description": "",
                "length_bp": 100,
                "sequence_sha256": haplotype_id,
            }
            for haplotype_id in haplotype_products
        ],
        "compiled_products": products,
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def write_molecules(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "read_id",
        "status",
        "assigned_product_id",
        "assignment_state",
        "best_sequence_sha256",
        "candidate_haplotypes",
        "edit_rate",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def molecule(read_id: str, product_id: str, sequence_hash: str, candidates: str) -> dict[str, str]:
    return {
        "read_id": read_id,
        "status": "complete",
        "assigned_product_id": product_id,
        "assignment_state": (
            "unique_sequence" if ";" not in candidates else "equivalent_haplotypes"
        ),
        "best_sequence_sha256": sequence_hash,
        "candidate_haplotypes": candidates,
        "edit_rate": "0.02",
    }


def run_model(compiled: Path, molecules: Path, config: GenotypeConfig | None = None):
    space = load_compiled_space(compiled)
    evidence = read_molecule_evidence(molecules)
    config = config or GenotypeConfig(
        min_assigned_reads=2,
        resolved_posterior=0.70,
        resolved_odds=2.0,
    )
    scores, effective = score_genotypes(space, evidence, config)
    return scores, make_call(scores, evidence, effective, config)


def test_resolves_pair_with_both_product_sequences(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(compiled, {"A": [("P", "ha")], "B": [("P", "hb")]})
    write_molecules(
        molecules,
        [
            *[molecule(f"a{i}", "P", "ha", "A") for i in range(5)],
            *[molecule(f"b{i}", "P", "hb", "B") for i in range(5)],
        ],
    )
    scores, call = run_model(compiled, molecules)
    assert call.call_state == "RESOLVED_RESEARCH_CALL"
    assert call.top_class is not None
    assert call.top_class.genotype_pairs == (("A", "B"),)
    assert scores[0].posterior > 0.70


def test_identical_haplotypes_are_reported_as_assay_equivalent(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(
        compiled,
        {"A": [("P", "ha")], "A_alias": [("P", "ha")], "B": [("Q", "hb")]},
    )
    write_molecules(
        molecules,
        [molecule(f"a{i}", "P", "ha", "A;A_alias") for i in range(8)],
    )
    _scores, call = run_model(compiled, molecules)
    assert call.call_state == "AMBIGUOUS_ASSAY_EQUIVALENT"
    assert call.top_class is not None
    assert set(call.top_class.genotype_pairs) == {
        ("A", "A"),
        ("A", "A_alias"),
        ("A_alias", "A_alias"),
    }


def test_one_allele_only_leaves_posterior_ambiguous(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(compiled, {"A": [("P", "ha")], "B": [("P", "hb")]})
    write_molecules(
        molecules,
        [molecule(f"a{i}", "P", "ha", "A") for i in range(4)],
    )
    _scores, call = run_model(compiled, molecules)
    assert call.call_state == "AMBIGUOUS_POSTERIOR"


def test_insufficient_evidence_is_no_call(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(compiled, {"A": [("P", "ha")], "B": [("P", "hb")]})
    write_molecules(molecules, [molecule("only", "P", "ha", "A")])
    _scores, call = run_model(compiled, molecules)
    assert call.call_state == "NO_CALL_INSUFFICIENT_EVIDENCE"


def test_synonymous_catalogue_labels_do_not_gain_prior_mass(tmp_path: Path) -> None:
    base_compiled = tmp_path / "base.json"
    expanded_compiled = tmp_path / "expanded.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(base_compiled, {"A": [("P", "ha")], "B": [("P", "hb")]})
    write_compiled(
        expanded_compiled,
        {"A": [("P", "ha")], "A_alias": [("P", "ha")], "B": [("P", "hb")]},
    )
    write_molecules(
        molecules,
        [
            *[molecule(f"a{i}", "P", "ha", "A;A_alias") for i in range(4)],
            *[molecule(f"b{i}", "P", "hb", "B") for i in range(4)],
        ],
    )
    base_space = load_compiled_space(base_compiled)
    expanded_space = load_compiled_space(expanded_compiled)
    evidence = read_molecule_evidence(molecules)
    config = GenotypeConfig(min_assigned_reads=2)
    base_scores, _ = score_genotypes(base_space, evidence, config)
    expanded_scores, _ = score_genotypes(expanded_space, evidence, config)

    base_ab = next(
        score for score in base_scores if score.genotype_pairs == (("A", "B"),)
    )
    expanded_ab = next(
        score
        for score in expanded_scores
        if set(score.genotype_pairs) == {("A", "B"), ("A_alias", "B")}
    )
    assert expanded_ab.posterior == pytest.approx(base_ab.posterior, rel=1e-12)


def test_cli_writes_call_and_ranked_posteriors(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(compiled, {"A": [("P", "ha")], "B": [("P", "hb")]})
    write_molecules(
        molecules,
        [
            *[molecule(f"a{i}", "P", "ha", "A") for i in range(5)],
            *[molecule(f"b{i}", "P", "hb", "B") for i in range(5)],
        ],
    )
    config = tmp_path / "config.yml"
    config.write_text(
        "min_assigned_reads: 2\nresolved_posterior: 0.70\nresolved_odds: 2.0\n",
        encoding="utf-8",
    )
    posteriors = tmp_path / "posteriors.tsv"
    call = tmp_path / "call.json"
    summary = tmp_path / "summary.json"
    assert genotype_main(
        [
            "--compiled-json",
            str(compiled),
            "--molecules-tsv",
            str(molecules),
            "--model-config",
            str(config),
            "--posteriors-tsv",
            str(posteriors),
            "--call-json",
            str(call),
            "--summary-json",
            str(summary),
        ]
    ) == 0
    assert json.loads(call.read_text(encoding="utf-8"))["call_state"] == "RESOLVED_RESEARCH_CALL"
    assert posteriors.read_text(encoding="utf-8").startswith("rank\tclass_id")
    assert (
        json.loads(summary.read_text(encoding="utf-8"))["model"]["calibration_status"]
        == "research_uncalibrated"
    )
