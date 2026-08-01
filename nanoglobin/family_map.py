"""HBA-family coordinate projection and copy-marker discovery.

The map is generated from sequence rather than maintained as a hand-written
population table. One sequence is selected as the family-coordinate anchor
(typically an HBA2-containing reference haplotype). HBA1, HBA2, hybrid and
structural-haplotype sequences are aligned to that anchor and projected into:

* reference-position keys (``R:<1-based-position>``); and
* insertion keys (``I:<number-of-anchor-bases-before-insertion>:<index>``).

This is sufficient to project copy-specific observations, identify paralogous
sequence variants, and expose intervals where an assay cannot distinguish copies.
Adding population haplotypes means adding sequences and regenerating the map; it
is not a separate manually curated "full population equivalence map" project.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from hashlib import sha256
from math import log2
from pathlib import Path
from typing import Any, Iterable, Mapping
import csv
import json
import re

from rapidfuzz.distance import Levenshtein


DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVN-]+$")


class FamilyMapError(ValueError):
    """Raised when family-coordinate inputs are invalid."""


@dataclass(frozen=True)
class SequenceRecord:
    sequence_id: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class CoordinateObservation:
    sequence_id: str
    family_key: str
    operation: str
    anchor_pos1: int | None
    insertion_after_pos1: int | None
    insertion_index: int | None
    sequence_pos1: int | None
    anchor_base: str
    sequence_base: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MarkerObservation:
    family_key: str
    marker_class: str
    anchor_pos1: int | None
    insertion_after_pos1: int | None
    insertion_index: int | None
    anchor_base: str
    allele_count: int
    entropy_bits: float
    observed_sequences: int
    missing_sequences: int
    unique_alleles: tuple[str, ...]
    alleles: Mapping[str, str]

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["unique_alleles"] = list(self.unique_alleles)
        value["alleles"] = dict(self.alleles)
        return value


def normalise_sequence(value: str, *, label: str = "sequence") -> str:
    sequence = re.sub(r"\s+", "", str(value)).upper().replace("U", "T")
    if not sequence:
        raise FamilyMapError(f"{label} is empty")
    if not DNA_RE.fullmatch(sequence):
        invalid = sorted(set(sequence) - set("ACGTRYSWKMBDHVN-"))
        raise FamilyMapError(
            f"{label} contains unsupported symbol(s): {', '.join(invalid)}"
        )
    if "-" in sequence:
        raise FamilyMapError(
            f"{label} contains alignment gaps; provide ungapped biological sequence"
        )
    return sequence


def read_fasta(path: str | Path) -> list[SequenceRecord]:
    """Read a FASTA catalogue with stable, unique identifiers."""

    records: list[SequenceRecord] = []
    current_id: str | None = None
    current_description = ""
    chunks: list[str] = []

    def flush() -> None:
        nonlocal current_id, current_description, chunks
        if current_id is None:
            return
        records.append(
            SequenceRecord(
                sequence_id=current_id,
                description=current_description,
                sequence=normalise_sequence(
                    "".join(chunks), label=f"FASTA record {current_id}"
                ),
            )
        )
        current_id = None
        current_description = ""
        chunks = []

    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                if not header:
                    raise FamilyMapError(f"empty FASTA header at line {line_number}")
                parts = header.split(maxsplit=1)
                current_id = parts[0]
                current_description = parts[1] if len(parts) == 2 else ""
            else:
                if current_id is None:
                    raise FamilyMapError(
                        f"sequence before first FASTA header at line {line_number}"
                    )
                chunks.append(line)
    flush()

    if not records:
        raise FamilyMapError(f"{path} contains no FASTA records")
    duplicates = sorted(
        sequence_id
        for sequence_id, count in Counter(r.sequence_id for r in records).items()
        if count > 1
    )
    if duplicates:
        raise FamilyMapError(f"duplicate FASTA identifiers: {', '.join(duplicates)}")
    return records


def _reference_key(pos0: int) -> str:
    return f"R:{pos0 + 1}"


def _insertion_key(anchor_bases_before: int, insertion_index: int) -> str:
    return f"I:{anchor_bases_before}:{insertion_index}"


def align_to_anchor(
    *,
    anchor: str,
    sequence: str,
    sequence_id: str,
) -> list[CoordinateObservation]:
    """Project one ungapped sequence into anchor and insertion coordinates."""

    anchor_n = normalise_sequence(anchor, label="anchor sequence")
    sequence_n = normalise_sequence(sequence, label=f"sequence {sequence_id}")
    editops = list(Levenshtein.editops(anchor_n, sequence_n))

    observations: list[CoordinateObservation] = []
    anchor_pos = 0
    sequence_pos = 0
    insertion_counts: Counter[int] = Counter()

    def emit_equal_until(target_anchor: int, target_sequence: int) -> None:
        nonlocal anchor_pos, sequence_pos
        while anchor_pos < target_anchor and sequence_pos < target_sequence:
            anchor_base = anchor_n[anchor_pos]
            source_base = sequence_n[sequence_pos]
            operation = "match" if anchor_base == source_base else "substitution"
            observations.append(
                CoordinateObservation(
                    sequence_id=sequence_id,
                    family_key=_reference_key(anchor_pos),
                    operation=operation,
                    anchor_pos1=anchor_pos + 1,
                    insertion_after_pos1=None,
                    insertion_index=None,
                    sequence_pos1=sequence_pos + 1,
                    anchor_base=anchor_base,
                    sequence_base=source_base,
                )
            )
            anchor_pos += 1
            sequence_pos += 1

    for edit in editops:
        emit_equal_until(edit.src_pos, edit.dest_pos)
        if edit.tag == "replace":
            observations.append(
                CoordinateObservation(
                    sequence_id=sequence_id,
                    family_key=_reference_key(anchor_pos),
                    operation="substitution",
                    anchor_pos1=anchor_pos + 1,
                    insertion_after_pos1=None,
                    insertion_index=None,
                    sequence_pos1=sequence_pos + 1,
                    anchor_base=anchor_n[anchor_pos],
                    sequence_base=sequence_n[sequence_pos],
                )
            )
            anchor_pos += 1
            sequence_pos += 1
        elif edit.tag == "delete":
            observations.append(
                CoordinateObservation(
                    sequence_id=sequence_id,
                    family_key=_reference_key(anchor_pos),
                    operation="deletion",
                    anchor_pos1=anchor_pos + 1,
                    insertion_after_pos1=None,
                    insertion_index=None,
                    sequence_pos1=None,
                    anchor_base=anchor_n[anchor_pos],
                    sequence_base="-",
                )
            )
            anchor_pos += 1
        elif edit.tag == "insert":
            insertion_counts[anchor_pos] += 1
            insertion_index = insertion_counts[anchor_pos]
            observations.append(
                CoordinateObservation(
                    sequence_id=sequence_id,
                    family_key=_insertion_key(anchor_pos, insertion_index),
                    operation="insertion",
                    anchor_pos1=None,
                    insertion_after_pos1=anchor_pos,
                    insertion_index=insertion_index,
                    sequence_pos1=sequence_pos + 1,
                    anchor_base="-",
                    sequence_base=sequence_n[sequence_pos],
                )
            )
            sequence_pos += 1
        else:  # pragma: no cover
            raise FamilyMapError(f"unsupported edit operation: {edit.tag}")

    emit_equal_until(len(anchor_n), len(sequence_n))
    while anchor_pos < len(anchor_n):
        observations.append(
            CoordinateObservation(
                sequence_id=sequence_id,
                family_key=_reference_key(anchor_pos),
                operation="deletion",
                anchor_pos1=anchor_pos + 1,
                insertion_after_pos1=None,
                insertion_index=None,
                sequence_pos1=None,
                anchor_base=anchor_n[anchor_pos],
                sequence_base="-",
            )
        )
        anchor_pos += 1
    while sequence_pos < len(sequence_n):
        insertion_counts[anchor_pos] += 1
        insertion_index = insertion_counts[anchor_pos]
        observations.append(
            CoordinateObservation(
                sequence_id=sequence_id,
                family_key=_insertion_key(anchor_pos, insertion_index),
                operation="insertion",
                anchor_pos1=None,
                insertion_after_pos1=anchor_pos,
                insertion_index=insertion_index,
                sequence_pos1=sequence_pos + 1,
                anchor_base="-",
                sequence_base=sequence_n[sequence_pos],
            )
        )
        sequence_pos += 1
    return observations


def build_coordinate_map(
    records: Iterable[SequenceRecord],
    *,
    anchor_id: str,
) -> tuple[SequenceRecord, list[CoordinateObservation]]:
    records_list = list(records)
    by_id = {record.sequence_id: record for record in records_list}
    if anchor_id not in by_id:
        raise FamilyMapError(f"anchor_id {anchor_id!r} is absent from the FASTA")
    if len(by_id) != len(records_list):
        raise FamilyMapError("sequence identifiers must be unique")
    anchor = by_id[anchor_id]

    observations: list[CoordinateObservation] = []
    for record in sorted(records_list, key=lambda item: item.sequence_id):
        observations.extend(
            align_to_anchor(
                anchor=anchor.sequence,
                sequence=record.sequence,
                sequence_id=record.sequence_id,
            )
        )
    return anchor, observations


def _entropy(values: Iterable[str]) -> float:
    counts = Counter(values)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((count / total) * log2(count / total) for count in counts.values())


def marker_rows(
    records: Iterable[SequenceRecord],
    observations: Iterable[CoordinateObservation],
    *,
    anchor_id: str,
) -> list[MarkerObservation]:
    """Summarize alleles and copy-informative markers at every family key."""

    records_list = sorted(records, key=lambda item: item.sequence_id)
    sequence_ids = [record.sequence_id for record in records_list]
    if anchor_id not in sequence_ids:
        raise FamilyMapError(f"anchor_id {anchor_id!r} is absent from records")
    by_key: dict[str, dict[str, CoordinateObservation]] = defaultdict(dict)
    for observation in observations:
        if observation.sequence_id in by_key[observation.family_key]:
            raise FamilyMapError(
                f"duplicate coordinate for {observation.sequence_id} at "
                f"{observation.family_key}"
            )
        by_key[observation.family_key][observation.sequence_id] = observation

    def sort_key(family_key: str) -> tuple[int, int, int]:
        kind, left, *rest = family_key.split(":")
        if kind == "R":
            return int(left), 0, 0
        return int(left), 1, int(rest[0])

    output: list[MarkerObservation] = []
    for family_key in sorted(by_key, key=sort_key):
        entries = by_key[family_key]
        exemplar = next(iter(entries.values()))
        alleles = {
            sequence_id: entries[sequence_id].sequence_base
            if sequence_id in entries
            else "-"
            for sequence_id in sequence_ids
        }
        concrete_values = [value for value in alleles.values() if value != "N"]
        allele_counts = Counter(concrete_values)
        allele_count = len(allele_counts)
        unique_alleles = tuple(
            sorted(
                f"{sequence_id}={allele}"
                for sequence_id, allele in alleles.items()
                if allele != "N" and allele_counts[allele] == 1
            )
        )
        if allele_count <= 1:
            marker_class = "invariant"
        elif "-" in allele_counts:
            marker_class = "presence_absence"
        elif unique_alleles:
            marker_class = "copy_informative"
        else:
            marker_class = "multi_allelic"
        output.append(
            MarkerObservation(
                family_key=family_key,
                marker_class=marker_class,
                anchor_pos1=exemplar.anchor_pos1,
                insertion_after_pos1=exemplar.insertion_after_pos1,
                insertion_index=exemplar.insertion_index,
                anchor_base=exemplar.anchor_base,
                allele_count=allele_count,
                entropy_bits=_entropy(concrete_values),
                observed_sequences=sum(value != "N" for value in alleles.values()),
                missing_sequences=sum(value == "N" for value in alleles.values()),
                unique_alleles=unique_alleles,
                alleles=alleles,
            )
        )
    return output


def summary_payload(
    *,
    anchor: SequenceRecord,
    records: Iterable[SequenceRecord],
    observations: Iterable[CoordinateObservation],
    markers: Iterable[MarkerObservation],
) -> dict[str, Any]:
    records_list = list(records)
    observations_list = list(observations)
    marker_list = list(markers)
    unique_by_sequence: Counter[str] = Counter()
    for marker in marker_list:
        for value in marker.unique_alleles:
            unique_by_sequence[value.split("=", 1)[0]] += 1
    return {
        "schema_version": 1,
        "anchor_id": anchor.sequence_id,
        "anchor_length_bp": len(anchor.sequence),
        "anchor_sha256": sha256(anchor.sequence.encode("ascii")).hexdigest(),
        "sequence_count": len(records_list),
        "sequence_ids": [record.sequence_id for record in records_list],
        "family_coordinate_count": len(marker_list),
        "reference_coordinate_count": sum(
            marker.family_key.startswith("R:") for marker in marker_list
        ),
        "insertion_coordinate_count": sum(
            marker.family_key.startswith("I:") for marker in marker_list
        ),
        "copy_informative_marker_count": sum(
            marker.marker_class
            in {"copy_informative", "presence_absence", "multi_allelic"}
            for marker in marker_list
        ),
        "unique_marker_count_by_sequence": dict(sorted(unique_by_sequence.items())),
        "operation_counts": dict(
            sorted(Counter(obs.operation for obs in observations_list).items())
        ),
    }


def write_positions_tsv(
    path: str | Path, observations: Iterable[CoordinateObservation]
) -> None:
    fields = [
        "sequence_id",
        "family_key",
        "operation",
        "anchor_pos1",
        "insertion_after_pos1",
        "insertion_index",
        "sequence_pos1",
        "anchor_base",
        "sequence_base",
    ]
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for observation in observations:
            writer.writerow(observation.to_dict())


def write_markers_tsv(path: str | Path, markers: Iterable[MarkerObservation]) -> None:
    fields = [
        "family_key",
        "marker_class",
        "anchor_pos1",
        "insertion_after_pos1",
        "insertion_index",
        "anchor_base",
        "allele_count",
        "entropy_bits",
        "observed_sequences",
        "missing_sequences",
        "unique_alleles",
        "alleles_json",
    ]
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for marker in markers:
            writer.writerow(
                {
                    "family_key": marker.family_key,
                    "marker_class": marker.marker_class,
                    "anchor_pos1": marker.anchor_pos1,
                    "insertion_after_pos1": marker.insertion_after_pos1,
                    "insertion_index": marker.insertion_index,
                    "anchor_base": marker.anchor_base,
                    "allele_count": marker.allele_count,
                    "entropy_bits": f"{marker.entropy_bits:.6f}",
                    "observed_sequences": marker.observed_sequences,
                    "missing_sequences": marker.missing_sequences,
                    "unique_alleles": ";".join(marker.unique_alleles),
                    "alleles_json": json.dumps(
                        dict(marker.alleles), sort_keys=True, separators=(",", ":")
                    ),
                }
            )


def write_map_json(
    path: str | Path,
    *,
    anchor: SequenceRecord,
    records: Iterable[SequenceRecord],
    observations: Iterable[CoordinateObservation],
    markers: Iterable[MarkerObservation],
) -> None:
    payload = {
        "schema_version": 1,
        "anchor": {
            "sequence_id": anchor.sequence_id,
            "description": anchor.description,
            "sequence": anchor.sequence,
            "sequence_sha256": sha256(anchor.sequence.encode("ascii")).hexdigest(),
        },
        "sequences": [asdict(record) for record in records],
        "coordinates": [observation.to_dict() for observation in observations],
        "markers": [marker.to_dict() for marker in markers],
    }
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


__all__ = [
    "CoordinateObservation",
    "FamilyMapError",
    "MarkerObservation",
    "SequenceRecord",
    "align_to_anchor",
    "build_coordinate_map",
    "marker_rows",
    "normalise_sequence",
    "read_fasta",
    "summary_payload",
    "write_map_json",
    "write_markers_tsv",
    "write_positions_tsv",
]
