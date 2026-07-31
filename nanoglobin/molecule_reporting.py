"""Streaming outputs and QC summaries for admitted amplicon molecules."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence
import csv

from .assay import AssayProfile
from .molecule_types import MoleculeObservation, OBSERVATION_FIELDS


PRODUCT_COUNT_FIELDS = [
    "product_id",
    "endpoint_candidate_reads",
    "assigned_reads",
    "unique_sequence_reads",
    "ambiguous_sequence_reads",
    "poor_fit_reads",
    "median_observed_length_bp",
    "min_observed_length_bp",
    "max_observed_length_bp",
    "fraction_of_all_reads",
    "fraction_of_complete_reads",
]


def write_observation_tsv_header(handle) -> csv.DictWriter:
    writer = csv.DictWriter(handle, fieldnames=OBSERVATION_FIELDS, delimiter="\t")
    writer.writeheader()
    return writer


def _weighted_median(counts: Counter[int]) -> float | None:
    """Return the exact median represented by an integer histogram."""

    total = sum(counts.values())
    if total == 0:
        return None
    lower_rank = (total - 1) // 2
    upper_rank = total // 2
    seen = 0
    lower_value: int | None = None
    upper_value: int | None = None
    for value, count in sorted(counts.items()):
        next_seen = seen + count
        if lower_value is None and lower_rank < next_seen:
            lower_value = value
        if upper_rank < next_seen:
            upper_value = value
            break
        seen = next_seen
    assert lower_value is not None and upper_value is not None
    return (lower_value + upper_value) / 2.0


class AdmissionAccumulator:
    """Streaming summary state with bounded integer length histograms."""

    def __init__(self) -> None:
        self.total = 0
        self.status_counts: Counter[str] = Counter()
        self.orientation_counts: Counter[str] = Counter()
        self.assignment_counts: Counter[str] = Counter()
        self.product_counts: Counter[str] = Counter()
        self.product_candidate_counts: Counter[str] = Counter()
        self.product_state_counts: dict[str, Counter[str]] = defaultdict(Counter)
        self.product_lengths: dict[str, Counter[int]] = defaultdict(Counter)

    def add(self, observation: MoleculeObservation) -> None:
        self.total += 1
        self.status_counts[observation.status] += 1
        self.orientation_counts[observation.orientation] += 1
        self.assignment_counts[observation.assignment_state or "none"] += 1
        for product_id in observation.candidate_product_ids:
            self.product_candidate_counts[product_id] += 1
        if observation.assigned_product_id:
            product_id = observation.assigned_product_id
            self.product_counts[product_id] += 1
            self.product_state_counts[product_id][
                observation.assignment_state or "none"
            ] += 1
            if observation.observed_product_length_bp is not None:
                self.product_lengths[product_id][
                    int(observation.observed_product_length_bp)
                ] += 1

    @property
    def complete_total(self) -> int:
        return self.status_counts["complete"]

    def product_rows(self) -> Iterator[dict[str, Any]]:
        product_ids = sorted(
            set(self.product_candidate_counts) | set(self.product_counts)
        )
        for product_id in product_ids:
            length_counts = self.product_lengths.get(product_id, Counter())
            median_length = _weighted_median(length_counts)
            states = self.product_state_counts.get(product_id, Counter())
            assigned = self.product_counts[product_id]
            yield {
                "product_id": product_id,
                "endpoint_candidate_reads": self.product_candidate_counts[product_id],
                "assigned_reads": assigned,
                "unique_sequence_reads": states["unique_sequence"],
                "ambiguous_sequence_reads": (
                    states["equivalent_haplotypes"] + states["low_margin"]
                ),
                "poor_fit_reads": states["poor_sequence_fit"],
                "median_observed_length_bp": (
                    f"{median_length:.1f}" if median_length is not None else ""
                ),
                "min_observed_length_bp": (
                    min(length_counts) if length_counts else ""
                ),
                "max_observed_length_bp": (
                    max(length_counts) if length_counts else ""
                ),
                "fraction_of_all_reads": (
                    f"{assigned / self.total:.6f}"
                    if self.total
                    else "0.000000"
                ),
                "fraction_of_complete_reads": (
                    f"{assigned / self.complete_total:.6f}"
                    if self.complete_total
                    else "0.000000"
                ),
            }

    def summary(self, *, sample: str, profile: AssayProfile) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "sample": sample,
            "assay_id": profile.assay_id,
            "assay_version": profile.version,
            "reads_total": self.total,
            "status_counts": dict(sorted(self.status_counts.items())),
            "orientation_counts": dict(sorted(self.orientation_counts.items())),
            "assignment_state_counts": dict(
                sorted(self.assignment_counts.items())
            ),
            "assigned_product_counts": dict(sorted(self.product_counts.items())),
            "complete_fraction": (
                self.complete_total / self.total if self.total else 0.0
            ),
            "off_target_fraction": (
                self.status_counts["off_target"] / self.total
                if self.total
                else 0.0
            ),
            "one_ended_fraction": (
                self.status_counts["one_ended"] / self.total
                if self.total
                else 0.0
            ),
            "chimera_candidate_fraction": (
                self.status_counts["chimera_candidate"] / self.total
                if self.total
                else 0.0
            ),
            "primer_dimer_candidate_fraction": (
                self.status_counts["primer_dimer_candidate"] / self.total
                if self.total
                else 0.0
            ),
            "unexpected_length_fraction": (
                self.status_counts["unexpected_length"] / self.total
                if self.total
                else 0.0
            ),
        }


def _as_accumulator(
    observations: Iterable[MoleculeObservation] | AdmissionAccumulator,
) -> AdmissionAccumulator:
    if isinstance(observations, AdmissionAccumulator):
        return observations
    accumulator = AdmissionAccumulator()
    for observation in observations:
        accumulator.add(observation)
    return accumulator


def write_product_counts_tsv(
    path: str | Path,
    observations: Iterable[MoleculeObservation] | AdmissionAccumulator,
) -> None:
    accumulator = _as_accumulator(observations)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=PRODUCT_COUNT_FIELDS, delimiter="\t"
        )
        writer.writeheader()
        for row in accumulator.product_rows():
            writer.writerow(row)


def summary_payload(
    *,
    sample: str,
    profile: AssayProfile,
    observations: Sequence[MoleculeObservation] | AdmissionAccumulator,
) -> dict[str, Any]:
    return _as_accumulator(observations).summary(sample=sample, profile=profile)
