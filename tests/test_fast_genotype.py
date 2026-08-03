from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from nanoglobin.fast_genotype import (
    cython_backend_available,
    score_genotypes_auto,
)
from nanoglobin.genotype import (
    GenotypeConfig,
    load_compiled_space,
    read_molecule_evidence,
)


def write_compiled(path: Path) -> None:
    haplotypes = {
        "A": [("P", "ha")],
        "B": [("P", "hb")],
        "C": [("Q", "hc")],
    }
    products = []
    for haplotype_id, keys in haplotypes.items():
        for product_id, sequence_hash in keys:
            products.append(
                {
                    "assay_id": "fast-test",
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
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "assay": {
                    "assay_id": "fast-test",
                    "version": "1",
                    "platform": "ont",
                },
                "haplotypes": [
                    {
                        "haplotype_id": haplotype_id,
                        "description": "",
                        "length_bp": 100,
                        "sequence_sha256": haplotype_id,
                    }
                    for haplotype_id in haplotypes
                ],
                "compiled_products": products,
            }
        ),
        encoding="utf-8",
    )


def write_molecules(path: Path) -> None:
    fields = [
        "read_id",
        "status",
        "assigned_product_id",
        "assignment_state",
        "best_sequence_sha256",
        "candidate_haplotypes",
        "edit_rate",
    ]
    rows = [
        {
            "read_id": "a_exact_1",
            "status": "complete",
            "assigned_product_id": "P",
            "assignment_state": "unique_sequence",
            "best_sequence_sha256": "ha",
            "candidate_haplotypes": "A",
            "edit_rate": "0.01",
        },
        {
            "read_id": "a_exact_2",
            "status": "complete",
            "assigned_product_id": "P",
            "assignment_state": "unique_sequence",
            "best_sequence_sha256": "ha",
            "candidate_haplotypes": "A",
            "edit_rate": "0.03",
        },
        {
            "read_id": "b_no_hash",
            "status": "complete",
            "assigned_product_id": "P",
            "assignment_state": "low_margin",
            "best_sequence_sha256": "",
            "candidate_haplotypes": "B",
            "edit_rate": "0.05",
        },
        {
            # A supplied but unmodelled sequence hash must receive the same
            # product-only quarter-probability as the Python reference, not the
            # no-hash/candidate-overlap half-probability.
            "read_id": "unknown_hash",
            "status": "complete",
            "assigned_product_id": "P",
            "assignment_state": "unique_sequence",
            "best_sequence_sha256": "not-in-catalogue",
            "candidate_haplotypes": "A",
            "edit_rate": "0.10",
        },
        {
            "read_id": "c_exact",
            "status": "complete",
            "assigned_product_id": "Q",
            "assignment_state": "unique_sequence",
            "best_sequence_sha256": "hc",
            "candidate_haplotypes": "C",
            "edit_rate": "0.02",
        },
        {
            "read_id": "off_target",
            "status": "off_target",
            "assigned_product_id": "",
            "assignment_state": "",
            "best_sequence_sha256": "",
            "candidate_haplotypes": "",
            "edit_rate": "",
        },
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_cython_backend_matches_python_reference(tmp_path: Path) -> None:
    if not cython_backend_available():
        pytest.skip("Cython extension was not built in this environment")
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(compiled)
    write_molecules(molecules)
    space = load_compiled_space(compiled)
    evidence = read_molecule_evidence(molecules)
    config = GenotypeConfig(
        min_assigned_reads=2,
        effective_read_cap=4.5,
        effective_count_cap=3.25,
        product_efficiencies={"P": 1.7, "Q": 0.6},
        dropout_probability=0.08,
        artifact_probability=0.03,
    )

    python_scores, python_effective, python_backend = score_genotypes_auto(
        space, evidence, config, backend="python"
    )
    cython_scores, cython_effective, cython_backend = score_genotypes_auto(
        space, evidence, config, backend="cython"
    )

    assert python_backend == "python_reference"
    assert cython_backend == "cython_dense_grouped"
    assert cython_effective == pytest.approx(python_effective, rel=1e-13, abs=1e-13)
    assert [score.class_id for score in cython_scores] == [
        score.class_id for score in python_scores
    ]
    for reference, accelerated in zip(python_scores, cython_scores):
        assert accelerated.genotype_pairs == reference.genotype_pairs
        assert accelerated.signature == reference.signature
        for field in (
            "log_read_likelihood",
            "log_count_likelihood",
            "log_dropout_likelihood",
            "log_prior",
            "log_score",
            "posterior",
        ):
            assert getattr(accelerated, field) == pytest.approx(
                getattr(reference, field), rel=1e-12, abs=1e-12
            )


def test_auto_backend_uses_compiled_kernel_when_available(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(compiled)
    write_molecules(molecules)
    space = load_compiled_space(compiled)
    evidence = read_molecule_evidence(molecules)
    _scores, _effective, backend = score_genotypes_auto(
        space, evidence, GenotypeConfig(min_assigned_reads=2), backend="auto"
    )
    if cython_backend_available():
        assert backend == "cython_dense_grouped"
    else:
        assert backend.startswith("python_reference")
