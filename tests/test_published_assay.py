from __future__ import annotations

import csv
from pathlib import Path

from nanoglobin.assay import load_assay_profile
from nanoglobin.calibration_io import read_control_manifest


ROOT = Path(__file__).resolve().parents[1]


def test_huang_2023_published_assay_profile_is_executable() -> None:
    profile = load_assay_profile(ROOT / "assays/huang_2023_ont_long_pcr.yml")
    assert profile.assay_id == "huang_2023_ont_long_pcr"
    assert profile.platform == "ont"
    assert len(profile.primers) == 12
    assert len(profile.products) == 8
    assert profile.primers["HBA_F"].sequence == "AGCTAGAGCATTGGTGGTCATGCCC"
    assert profile.primers["HBA_R"].sequence == "GACTTCGCGGTGGCTCCACTTTCC"
    assert profile.primers["HBD_F"].sequence == "ATGCTCATGGGTGTGTATTTGTCTGCC"
    assert profile.primers["HBB_R"].sequence == "ACAGAGACAACTAAGGCTGAGTGGCC"
    products = {item.product_id: item for item in profile.products}
    assert products["HBA_whole_and_rearranged"].forward_primer == "HBA_F"
    assert products["HBA_whole_and_rearranged"].reverse_primer == "HBA_R"
    assert products["alpha_27_6_deletion"].forward_primer == "TF_F"
    assert products["alpha_27_6_deletion"].reverse_primer == "HBA_R"
    assert products["SEA_HPFH_deletion"].forward_primer == "HBD_F"
    assert products["SEA_HPFH_deletion"].reverse_primer == "SH_R"
    assert products["Chinese_Ggamma_Agammadeltabeta0_deletion"].forward_primer == (
        "HBG_F"
    )
    assert products["Chinese_Ggamma_Agammadeltabeta0_deletion"].reverse_primer == (
        "DB_R"
    )
    assert all(item.required for item in profile.products)


def test_validation_resource_registry_is_complete_and_unique() -> None:
    path = ROOT / "validation/resources.tsv"
    required = {
        "resource_id",
        "resource_type",
        "access_status",
        "assay_match",
        "truth_role",
        "genotype_scope",
        "identifier",
        "intended_use",
        "next_action",
        "critical_caveat",
        "source",
    }
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert set(reader.fieldnames or []) == required
        rows = list(reader)
    assert len(rows) >= 10
    resource_ids = [row["resource_id"] for row in rows]
    assert len(resource_ids) == len(set(resource_ids))
    assert all(all((row[field] or "").strip() for field in required) for row in rows)
    by_id = {row["resource_id"]: row for row in rows}
    assert by_id["PRJNA1439314"]["truth_role"] == "no_independent_truth"
    assert by_id["DRAGEN_selected_HBA"]["truth_role"] == "comparator_only"
    assert by_id["internal_amplidex_paired"]["assay_match"] == (
        "exact_target_assay"
    )


def test_control_manifest_template_uses_declared_truth_scopes() -> None:
    rows = read_control_manifest(ROOT / "validation/control_manifest.template.tsv")
    assert len(rows) == 3
    assert {item.truth_status for item in rows} == {
        "orthogonally_resolved",
        "certified_reference",
        "synthetic",
    }
    assert {item.assay_scope for item in rows} == {
        "assay_matched",
        "synthetic_assay",
    }
