from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from nanoglobin.calibration_evidence import build_control_evidence
from nanoglobin.calibration_fit import calibrate_model
from nanoglobin.calibration_io import (
    load_calibration_provenance,
    read_control_manifest,
)
from nanoglobin.calibration_types import CalibrationError, CalibrationSettings
from nanoglobin.genotype import GenotypeConfig, load_compiled_space, load_config
from scripts.calibrate_amplicon_model import main as calibrate_main


def write_compiled(
    path: Path, haplotype_products: dict[str, list[tuple[str, str, bool]]]
) -> None:
    products = []
    for haplotype_id, keys in haplotype_products.items():
        for product_id, sequence_hash, required in keys:
            products.append(
                {
                    "assay_id": "cal-test",
                    "assay_version": "1",
                    "product_id": product_id,
                    "pool": "globin",
                    "roles": ["whole_gene"],
                    "required": required,
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
                    "assay_id": "cal-test",
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
                    for haplotype_id in haplotype_products
                ],
                "compiled_products": products,
            }
        ),
        encoding="utf-8",
    )


def write_molecules(
    path: Path,
    counts: dict[tuple[str, str, str], int],
    *,
    unsupported: int = 0,
) -> None:
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
    for (product_id, sequence_hash, candidates), count in counts.items():
        for index in range(count):
            rows.append(
                {
                    "read_id": f"{product_id}_{sequence_hash}_{index}",
                    "status": "complete",
                    "assigned_product_id": product_id,
                    "assignment_state": "unique_sequence",
                    "best_sequence_sha256": sequence_hash,
                    "candidate_haplotypes": candidates,
                    "edit_rate": "0.01",
                }
            )
    for index in range(unsupported):
        rows.append(
            {
                "read_id": f"unsupported_{index}",
                "status": "complete",
                "assigned_product_id": "P",
                "assignment_state": "unique_sequence",
                "best_sequence_sha256": "unexpected",
                "candidate_haplotypes": "A",
                "edit_rate": "0.10",
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "sample_id",
        "molecules_tsv",
        "haplotype_1",
        "haplotype_2",
        "stratum",
        "truth_status",
        "truth_source",
        "assay_scope",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_recovers_relative_product_efficiency(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    write_compiled(
        compiled,
        {"A": [("P", "hp", True), ("Q", "hq", True)]},
    )
    controls = []
    manifest_rows = []
    for sample_index in range(3):
        molecules = tmp_path / f"control_{sample_index}.tsv"
        write_molecules(
            molecules,
            {
                ("P", "hp", "A"): 80,
                ("Q", "hq", "A"): 20,
            },
        )
        controls.append(molecules)
        manifest_rows.append(
            {
                "sample_id": f"C{sample_index}",
                "molecules_tsv": molecules.name,
                "haplotype_1": "A",
                "haplotype_2": "A",
                "stratum": "lot-1",
                "truth_status": "orthogonally_resolved",
                "truth_source": "gap-PCR plus junction sequencing",
                "assay_scope": "assay_matched",
            }
        )
    manifest = tmp_path / "controls.tsv"
    write_manifest(manifest, manifest_rows)

    space = load_compiled_space(compiled)
    rows = read_control_manifest(manifest)
    result = calibrate_model(
        space,
        rows,
        base_config=GenotypeConfig(
            effective_count_cap=1000.0,
            effective_read_cap=1000.0,
        ),
        settings=CalibrationSettings(
            min_effective_reads=1.0,
            efficiency_prior_strength=1e-6,
            min_controls_per_stratum=2,
        ),
    )

    efficiencies = result.calibrated_config.product_efficiencies
    assert efficiencies["P"] / efficiencies["Q"] == pytest.approx(4.0, rel=1e-5)
    assert result.calibration_status == "assay_matched_control_calibrated_research"
    assert result.global_dropout["dropout_events"] == 0
    assert result.global_dropout["dropout_trials"] == 6
    assert result.strata["lot-1"]["status"] == "empirical_bayes_stratum_fit"


def test_comparator_only_controls_are_excluded_by_default(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "control.tsv"
    manifest = tmp_path / "controls.tsv"
    write_compiled(compiled, {"A": [("P", "hp", True)]})
    write_molecules(molecules, {("P", "hp", "A"): 20})
    write_manifest(
        manifest,
        [
            {
                "sample_id": "DRAGEN_SELECTED",
                "molecules_tsv": molecules.name,
                "haplotype_1": "A",
                "haplotype_2": "A",
                "stratum": "global",
                "truth_status": "comparator_only",
                "truth_source": "DRAGEN HBA caller",
                "assay_scope": "caller_comparator",
            }
        ],
    )
    space = load_compiled_space(compiled)
    rows = read_control_manifest(manifest)
    evidence = build_control_evidence(
        space,
        rows,
        GenotypeConfig(),
        CalibrationSettings(min_effective_reads=1.0),
    )
    assert len(evidence) == 1
    assert not evidence[0].eligible
    assert "comparator_only" in evidence[0].exclusion_reason
    with pytest.raises(CalibrationError, match="no eligible calibration controls"):
        calibrate_model(
            space,
            rows,
            settings=CalibrationSettings(min_effective_reads=1.0),
        )


def test_disconnected_single_product_components_remain_unidentified(
    tmp_path: Path,
) -> None:
    compiled = tmp_path / "compiled.json"
    a_reads = tmp_path / "a.tsv"
    b_reads = tmp_path / "b.tsv"
    manifest = tmp_path / "controls.tsv"
    write_compiled(
        compiled,
        {
            "A": [("P", "hp", True)],
            "B": [("Q", "hq", True)],
        },
    )
    write_molecules(a_reads, {("P", "hp", "A"): 50})
    write_molecules(b_reads, {("Q", "hq", "B"): 10})
    write_manifest(
        manifest,
        [
            {
                "sample_id": "A_CONTROL",
                "molecules_tsv": a_reads.name,
                "haplotype_1": "A",
                "haplotype_2": "A",
                "stratum": "global",
                "truth_status": "synthetic",
                "truth_source": "sealed fixture A",
                "assay_scope": "synthetic_assay",
            },
            {
                "sample_id": "B_CONTROL",
                "molecules_tsv": b_reads.name,
                "haplotype_1": "B",
                "haplotype_2": "B",
                "stratum": "global",
                "truth_status": "synthetic",
                "truth_source": "sealed fixture B",
                "assay_scope": "synthetic_assay",
            },
        ],
    )
    result = calibrate_model(
        load_compiled_space(compiled),
        read_control_manifest(manifest),
        settings=CalibrationSettings(min_effective_reads=1.0),
    )
    assert len(result.components) == 2
    assert all(
        item.identifiability == "single_product_relative_efficiency_unidentifiable"
        for item in result.components
    )
    assert any("scale is not identified" in warning for warning in result.warnings)
    assert result.calibration_status == "synthetic_fitted_research"


def test_calibration_cli_writes_provenance_checked_config(tmp_path: Path) -> None:
    compiled = tmp_path / "compiled.json"
    molecules = tmp_path / "control.tsv"
    manifest = tmp_path / "controls.tsv"
    write_compiled(
        compiled,
        {"A": [("P", "hp", True), ("Q", "hq", True)]},
    )
    write_molecules(
        molecules,
        {("P", "hp", "A"): 30, ("Q", "hq", "A"): 15},
        unsupported=2,
    )
    write_manifest(
        manifest,
        [
            {
                "sample_id": "CONTROL",
                "molecules_tsv": molecules.name,
                "haplotype_1": "A",
                "haplotype_2": "A",
                "stratum": "lot-1",
                "truth_status": "certified_reference",
                "truth_source": "reference material certificate",
                "assay_scope": "assay_matched",
            }
        ],
    )
    calibration_json = tmp_path / "calibration.json"
    model_config = tmp_path / "calibrated.yml"
    product_tsv = tmp_path / "products.tsv"
    fit_tsv = tmp_path / "fit.tsv"
    assert calibrate_main(
        [
            "--compiled-json",
            str(compiled),
            "--control-manifest",
            str(manifest),
            "--calibration-json",
            str(calibration_json),
            "--calibrated-model-config",
            str(model_config),
            "--product-parameters-tsv",
            str(product_tsv),
            "--sample-product-fit-tsv",
            str(fit_tsv),
        ]
    ) == 0
    space = load_compiled_space(compiled)
    config = load_config(model_config)
    provenance = load_calibration_provenance(calibration_json, space, config)
    assert provenance is not None
    assert provenance["calibration_status"] == (
        "assay_matched_control_calibrated_research"
    )
    assert product_tsv.read_text(encoding="utf-8").startswith(
        "component_id\tproduct_id"
    )
    assert "CONTROL" in fit_tsv.read_text(encoding="utf-8")

    values = json.loads(calibration_json.read_text(encoding="utf-8"))
    values["calibrated_model_config"]["dropout_probability"] = 0.999
    calibration_json.write_text(json.dumps(values), encoding="utf-8")
    with pytest.raises(CalibrationError, match="does not match"):
        load_calibration_provenance(calibration_json, space, config)
