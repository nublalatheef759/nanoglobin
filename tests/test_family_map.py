from __future__ import annotations

import json
from pathlib import Path

from nanoglobin.family_map import (
    SequenceRecord,
    align_to_anchor,
    build_coordinate_map,
    marker_rows,
    summary_payload,
)
from scripts.build_hba_family_map import main as family_map_main


def test_alignment_projects_substitution_deletion_and_insertion() -> None:
    anchor = "ACGTACGT"
    substitution = align_to_anchor(
        anchor=anchor, sequence="ACCTACGT", sequence_id="HBA1"
    )
    assert any(
        row.family_key == "R:3"
        and row.operation == "substitution"
        and row.anchor_base == "G"
        and row.sequence_base == "C"
        for row in substitution
    )

    deletion = align_to_anchor(anchor=anchor, sequence="ACGACGT", sequence_id="DEL")
    assert any(
        row.operation == "deletion" and row.anchor_base == "T"
        for row in deletion
    )

    insertion = align_to_anchor(
        anchor=anchor, sequence="ACGTTACGT", sequence_id="INS"
    )
    assert any(
        row.family_key.startswith("I:4:")
        and row.operation == "insertion"
        and row.sequence_base == "T"
        for row in insertion
    )


def test_marker_map_is_generated_from_added_sequences() -> None:
    records = [
        SequenceRecord("HBA2", "ACGTACGT"),
        SequenceRecord("HBA1", "ACCTACGT"),
        SequenceRecord("HYBRID", "ACGACGT"),
    ]
    anchor, positions = build_coordinate_map(records, anchor_id="HBA2")
    markers = marker_rows(records, positions, anchor_id="HBA2")

    informative = {
        marker.family_key: marker
        for marker in markers
        if marker.marker_class != "invariant"
    }
    assert informative
    assert any(
        marker.marker_class in {"copy_informative", "presence_absence"}
        for marker in informative.values()
    )

    summary = summary_payload(
        anchor=anchor,
        records=records,
        observations=positions,
        markers=markers,
    )
    assert summary["sequence_count"] == 3
    assert summary["copy_informative_marker_count"] >= 1


def test_family_map_cli_writes_inspectable_outputs(tmp_path: Path) -> None:
    fasta = tmp_path / "copies.fa"
    fasta.write_text(
        ">HBA2 anchor\nACGTACGT\n>HBA1 copy\nACCTACGT\n>HYBRID\nACGACGT\n",
        encoding="utf-8",
    )
    positions = tmp_path / "positions.tsv"
    markers = tmp_path / "markers.tsv"
    map_json = tmp_path / "map.json"
    summary = tmp_path / "summary.json"

    assert family_map_main(
        [
            "--sequences",
            str(fasta),
            "--anchor-id",
            "HBA2",
            "--positions-tsv",
            str(positions),
            "--markers-tsv",
            str(markers),
            "--map-json",
            str(map_json),
            "--summary-json",
            str(summary),
        ]
    ) == 0
    assert positions.exists() and markers.exists() and map_json.exists()
    payload = json.loads(summary.read_text(encoding="utf-8"))
    assert payload["anchor_id"] == "HBA2"
    assert payload["sequence_count"] == 3
