"""Typed records for primer-aware amplicon molecule evidence."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from .assay import CompiledProduct


class MoleculeAdmissionError(ValueError):
    """Raised when FASTQ, compiled-assay or admission input is inconsistent."""


@dataclass(frozen=True)
class FastqRecord:
    read_id: str
    sequence: str
    qualities: str


@dataclass(frozen=True)
class AdmissionConfig:
    """Read-evidence thresholds, distinct from genomic primer binding rules."""

    end_search_bp: int = 100
    primer_seed_length: int = 7
    min_primer_overlap: int = 10
    max_primer_mismatches: int = 4
    max_primer_error_rate: float = 0.25
    terminal_offset_penalty: float = 0.08
    clipped_primer_penalty: float = 1.25
    length_tolerance_fraction: float = 0.10
    max_sequence_edit_rate: float = 0.25
    sequence_assignment_margin: float = 0.02
    assumed_sequence_error_rate: float = 0.10

    def validate(self) -> None:
        if self.end_search_bp < 1:
            raise MoleculeAdmissionError("end_search_bp must be >= 1")
        if self.primer_seed_length < 4:
            raise MoleculeAdmissionError("primer_seed_length must be >= 4")
        if self.min_primer_overlap < 4:
            raise MoleculeAdmissionError("min_primer_overlap must be >= 4")
        if self.max_primer_mismatches < 0:
            raise MoleculeAdmissionError("max_primer_mismatches must be >= 0")
        if not 0 <= self.max_primer_error_rate <= 1:
            raise MoleculeAdmissionError(
                "max_primer_error_rate must be between 0 and 1"
            )
        if self.terminal_offset_penalty < 0 or self.clipped_primer_penalty < 0:
            raise MoleculeAdmissionError("primer penalties must be >= 0")
        if not 0 <= self.length_tolerance_fraction <= 1:
            raise MoleculeAdmissionError(
                "length_tolerance_fraction must be between 0 and 1"
            )
        if not 0 <= self.max_sequence_edit_rate <= 1:
            raise MoleculeAdmissionError(
                "max_sequence_edit_rate must be between 0 and 1"
            )
        if self.sequence_assignment_margin < 0:
            raise MoleculeAdmissionError("sequence_assignment_margin must be >= 0")
        if not 0 < self.assumed_sequence_error_rate < 0.75:
            raise MoleculeAdmissionError(
                "assumed_sequence_error_rate must be in (0, 0.75)"
            )


@dataclass(frozen=True)
class TerminalPrimerMatch:
    primer_id: str
    side: str
    read_start0: int
    read_end0: int
    primer_start0: int
    primer_end0: int
    primer_length_bp: int
    overlap_bp: int
    mismatches: int
    offset_bp: int
    quality_log_likelihood: float
    score: float

    @property
    def clipped_bases(self) -> int:
        return self.primer_start0 + (self.primer_length_bp - self.primer_end0)


@dataclass(frozen=True)
class SequenceCompatibility:
    product_id: str
    sequence_sha256: str
    haplotype_ids: tuple[str, ...]
    edit_distance: int
    edit_rate: float
    log_likelihood: float
    target_length_bp: int
    observed_length_bp: int


@dataclass(frozen=True)
class MoleculeObservation:
    sample: str
    read_id: str
    read_length_bp: int
    status: str
    orientation: str
    left_primer_id: str = ""
    right_primer_id: str = ""
    left_overlap_bp: int | None = None
    right_overlap_bp: int | None = None
    left_mismatches: int | None = None
    right_mismatches: int | None = None
    trim_start0: int | None = None
    trim_end0: int | None = None
    observed_product_length_bp: int | None = None
    candidate_product_ids: tuple[str, ...] = ()
    assigned_product_id: str = ""
    assignment_state: str = ""
    best_sequence_sha256: str = ""
    candidate_haplotypes: tuple[str, ...] = ()
    edit_distance: int | None = None
    edit_rate: float | None = None
    sequence_log_likelihood: float | None = None
    likelihood_margin: float | None = None
    reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["candidate_product_ids"] = ";".join(self.candidate_product_ids)
        value["candidate_haplotypes"] = ";".join(self.candidate_haplotypes)
        return value


@dataclass(frozen=True)
class CompiledAssay:
    assay_id: str
    assay_version: str
    products: tuple[CompiledProduct, ...]


OBSERVATION_FIELDS = [
    "sample",
    "read_id",
    "read_length_bp",
    "status",
    "orientation",
    "left_primer_id",
    "right_primer_id",
    "left_overlap_bp",
    "right_overlap_bp",
    "left_mismatches",
    "right_mismatches",
    "trim_start0",
    "trim_end0",
    "observed_product_length_bp",
    "candidate_product_ids",
    "assigned_product_id",
    "assignment_state",
    "best_sequence_sha256",
    "candidate_haplotypes",
    "edit_distance",
    "edit_rate",
    "sequence_log_likelihood",
    "likelihood_margin",
    "reason",
]
