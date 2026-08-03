from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from nanoglobin.genotype import GenotypeConfig, load_compiled_space, read_molecule_evidence
from nanoglobin.stability import (
    StabilityError,
    StabilityGrid,
    load_stability_grid,
    run_stability_grid,
    stability_summary,
)
from scripts.sweep_genotype_parameters import main as sweep_main


def write_compiled(path: Path) -> None:
    products = []
    for haplotype_id, sequence_hash in (("A", "ha"), ("B", "hb")):
        products.append(
            {
                "assay_id": "stability-test",
                "assay_version": "1",
                "product_id": "P",
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
                    "assay_id": "stability-test",
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
                    for haplotype_id in ("A", "B")
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
    rows = []
    for haplotype_id, sequence_hash in (("A", "ha"), ("B", "hb")):
        for index in range(8):
            rows.append(
                {
                    "read_id": f"{haplotype_id}_{index}",
                    "status": "complete",
                    "assigned_product_id": "P",
                    "assignment_state": "unique_sequence",
                    "best_sequence_sha256": sequence_hash,
                    "candidate_haplotypes": haplotype_id,
                    "edit_rate": "0.02",
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_stability_grid_reports_consistent_top_class(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    write_compiled(compiled)
    write_molecules(molecules)
    space = load_compiled_space(compiled)
    evidence = read_molecule_evidence(molecules)
    config = GenotypeConfig(
        min_assigned_reads=2,
        resolved_posterior=0.70,
        resolved_odds=2.0,
    )
    grid = StabilityGrid(
        axes={
            "dropout_probability": (0.01, 0.05, 0.20),
            "artifact_probability": (0.01, 0.05),
            "count_concentration": (5.0, 20.0),
        }
    )
    scenarios = run_stability_grid(
        space,
        evidence,
        config,
        grid,
    )
    summary = stability_summary(scenarios, base_config=config, grid=grid)
    assert len(scenarios) == 12
    assert summary["scenario_count"] == 12
    assert summary["modal_top_class_fraction"] == 1.0
    assert set(summary["call_state_counts"]) == {"RESOLVED_RESEARCH_CALL"}


def test_invalid_grid_is_rejected() -> None:
    config = GenotypeConfig()
    with pytest.raises(StabilityError, match="mapping-valued axes"):
        StabilityGrid(axes={"product_efficiencies": (1.0, 2.0)}).validate(config)
    with pytest.raises(StabilityError, match="maximum"):
        StabilityGrid(
            axes={
                "dropout_probability": tuple(
                    (index + 1) / 100 for index in range(20)
                ),
                "artifact_probability": tuple(
                    (index + 1) / 1000 for index in range(20)
                ),
            },
            max_scenarios=100,
        ).validate(config)


def test_stability_cli_writes_scenarios_and_summary(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "molecules.tsv"
    config = tmp_path / "config.yml"
    grid = tmp_path / "grid.yml"
    scenarios = tmp_path / "scenarios.tsv"
    summary = tmp_path / "summary.json"
    write_compiled(compiled)
    write_molecules(molecules)
    config.write_text(
        "min_assigned_reads: 2\nresolved_posterior: 0.70\nresolved_odds: 2.0\n",
        encoding="utf-8",
    )
    grid.write_text(
        "axes:\n"
        "  dropout_probability: [0.01, 0.05]\n"
        "  artifact_probability: [0.01, 0.03]\n",
        encoding="utf-8",
    )
    assert sweep_main(
        [
            "--compiled-json",
            str(compiled),
            "--molecules-tsv",
            str(molecules),
            "--model-config",
            str(config),
            "--grid",
            str(grid),
            "--scenarios-tsv",
            str(scenarios),
            "--summary-json",
            str(summary),
        ]
    ) == 0
    assert scenarios.read_text(encoding="utf-8").startswith(
        "scenario_id\tartifact_probability\tdropout_probability"
    )
    payload = json.loads(summary.read_text(encoding="utf-8"))
    assert payload["scenario_count"] == 4
    assert payload["modal_top_class_fraction"] == 1.0


def test_load_grid_preserves_integer_axes(tmp_path: Path) -> None:
    profile = tmp_path / "grid.yml"
    profile.write_text(
        "axes:\n  min_assigned_reads: [1, 2, 3]\n",
        encoding="utf-8",
    )
    grid = load_stability_grid(profile, GenotypeConfig())
    assert grid.axes["min_assigned_reads"] == (1, 2, 3)
