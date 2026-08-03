from __future__ import annotations

import csv
import json
from math import log
from pathlib import Path

import pytest

from nanoglobin.genotype import (
    GenotypeConfig,
    load_compiled_space,
    read_molecule_evidence,
    score_genotypes,
)


def write_single_product_space(path: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "assay": {
                    "assay_id": "likelihood-semantics",
                    "version": "1",
                    "platform": "ont",
                },
                "haplotypes": [
                    {
                        "haplotype_id": "A",
                        "description": "",
                        "length_bp": 100,
                        "sequence_sha256": "catalogue-A",
                    }
                ],
                "compiled_products": [
                    {
                        "assay_id": "likelihood-semantics",
                        "assay_version": "1",
                        "product_id": "P",
                        "pool": "alpha",
                        "roles": ["whole_gene"],
                        "required": True,
                        "haplotype_id": "A",
                        "start0": 0,
                        "end0": 100,
                        "length_bp": 100,
                        "forward_mismatches": 0,
                        "reverse_mismatches": 0,
                        "sequence": "A" * 100,
                        "sequence_sha256": "compiled-A",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )


def write_molecule(path: Path, sequence_hash: str) -> None:
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
        writer.writerow(
            {
                "read_id": "read-1",
                "status": "complete",
                "assigned_product_id": "P",
                "assignment_state": "unique_sequence",
                "best_sequence_sha256": sequence_hash,
                "candidate_haplotypes": "A",
                "edit_rate": "0.02",
            }
        )


def test_supplied_unmodelled_hash_uses_product_only_quarter_probability(
    tmp_path: Path,
) -> None:
    compiled = tmp_path / "compiled.json"
    unmodelled = tmp_path / "unmodelled.tsv"
    no_hash = tmp_path / "no-hash.tsv"
    write_single_product_space(compiled)
    write_molecule(unmodelled, "not-in-candidate-catalogue")
    write_molecule(no_hash, "")

    space = load_compiled_space(compiled)
    config = GenotypeConfig(
        min_assigned_reads=1,
        artifact_probability=0.02,
        effective_read_cap=10.0,
        effective_count_cap=10.0,
    )
    unmodelled_scores, _ = score_genotypes(
        space,
        read_molecule_evidence(unmodelled),
        config,
    )
    no_hash_scores, _ = score_genotypes(
        space,
        read_molecule_evidence(no_hash),
        config,
    )

    quarter_probability = max(
        config.artifact_probability,
        0.25 * (1.0 - config.artifact_probability),
    )
    half_probability = max(
        config.artifact_probability,
        0.50 * (1.0 - config.artifact_probability),
    )
    edit_penalty = 0.02

    assert unmodelled_scores[0].log_read_likelihood == pytest.approx(
        log(quarter_probability) - edit_penalty,
        rel=1e-12,
        abs=1e-12,
    )
    assert no_hash_scores[0].log_read_likelihood == pytest.approx(
        log(half_probability) - edit_penalty,
        rel=1e-12,
        abs=1e-12,
    )
    assert (
        no_hash_scores[0].log_read_likelihood
        - unmodelled_scores[0].log_read_likelihood
    ) == pytest.approx(
        log(half_probability / quarter_probability),
        rel=1e-12,
        abs=1e-12,
    )
