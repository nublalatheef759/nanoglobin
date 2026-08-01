"""Primer-aware ONT amplicon molecule admission and product evidence.

The assay compiler answers which products a chromosome haplotype can generate.
This module addresses the next physical layer: which declared product, if any,
can explain each observed FASTQ molecule?

The implementation deliberately stops before a diploid genotype call. It emits
inspectable molecule-level evidence: terminal primer matches, molecule state,
product assignment, sequence compatibility, ambiguity and failure reasons.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, Iterable, Mapping, Sequence

from .assay import (
    AssayProfile,
    CompiledProduct,
    ProductDefinition,
    normalise_dna,
    reverse_complement,
)
from .fastq import read_fastq
from .molecule_reporting import (
    AdmissionAccumulator,
    summary_payload,
    write_observation_tsv_header,
    write_product_counts_tsv,
)
from .molecule_types import (
    AdmissionConfig,
    CompiledAssay,
    FastqRecord,
    MoleculeAdmissionError,
    MoleculeObservation,
    OBSERVATION_FIELDS,
    TerminalPrimerMatch,
)
from .primer_matching import best_terminal_primer_match
from .product_evidence import (
    CompiledSequenceIndex,
    build_compiled_sequence_index,
    load_compiled_assay,
    score_compiled_sequences,
)


@dataclass(frozen=True)
class _EndpointCandidate:
    product: ProductDefinition
    orientation: str
    left: TerminalPrimerMatch
    right: TerminalPrimerMatch
    endpoint_score: float
    observed_start0: int
    observed_end0: int
    observed_length_bp: int


def _reverse_qualities(qualities: str) -> str:
    return qualities[::-1]


def _endpoint_candidates_for_orientation(
    profile: AssayProfile,
    sequence: str,
    qualities: str,
    orientation: str,
    config: AdmissionConfig,
    *,
    left_primer_ids: Iterable[str] | None = None,
    right_primer_ids: Iterable[str] | None = None,
) -> tuple[
    list[_EndpointCandidate],
    dict[str, TerminalPrimerMatch],
    dict[str, TerminalPrimerMatch],
]:
    """Match terminal primers and construct length-compatible product candidates."""

    left_matches: dict[str, TerminalPrimerMatch] = {}
    right_matches: dict[str, TerminalPrimerMatch] = {}
    left_ids = (
        tuple(left_primer_ids)
        if left_primer_ids is not None
        else tuple(profile.primers)
    )
    right_ids = (
        tuple(right_primer_ids)
        if right_primer_ids is not None
        else tuple(profile.primers)
    )

    for primer_id in left_ids:
        primer = profile.primers[primer_id]
        left = best_terminal_primer_match(
            sequence,
            qualities,
            primer_id=primer_id,
            primer_sequence=primer.sequence,
            side="left",
            config=config,
        )
        if left is not None:
            left_matches[primer_id] = left

    for primer_id in right_ids:
        primer = profile.primers[primer_id]
        right = best_terminal_primer_match(
            sequence,
            qualities,
            primer_id=primer_id,
            primer_sequence=reverse_complement(primer.sequence),
            side="right",
            config=config,
        )
        if right is not None:
            right_matches[primer_id] = right

    candidates: list[_EndpointCandidate] = []
    for product in profile.products:
        left = left_matches.get(product.forward_primer)
        right = right_matches.get(product.reverse_primer)
        if left is None or right is None or left.read_end0 > right.read_start0:
            continue

        observed_start = left.read_start0
        observed_end = right.read_end0
        observed_length = observed_end - observed_start
        left_missing = left.primer_start0
        reverse_pattern_length = len(
            profile.primers[product.reverse_primer].sequence
        )
        right_missing = reverse_pattern_length - right.primer_end0
        expected_min = max(1, product.min_length - left_missing - right_missing)
        expected_max = max(
            expected_min,
            product.max_length - left_missing - right_missing,
        )
        lower_tolerance = max(
            2, int(expected_min * config.length_tolerance_fraction)
        )
        upper_tolerance = max(
            2, int(expected_max * config.length_tolerance_fraction)
        )
        if not (
            expected_min - lower_tolerance
            <= observed_length
            <= expected_max + upper_tolerance
        ):
            continue

        midpoint = (expected_min + expected_max) / 2.0
        length_penalty = abs(observed_length - midpoint) / max(1.0, midpoint)
        candidates.append(
            _EndpointCandidate(
                product=product,
                orientation=orientation,
                left=left,
                right=right,
                endpoint_score=left.score + right.score - length_penalty,
                observed_start0=observed_start,
                observed_end0=observed_end,
                observed_length_bp=observed_length,
            )
        )

    candidates.sort(
        key=lambda item: (-item.endpoint_score, item.product.product_id)
    )
    return candidates, left_matches, right_matches


def _best_declared_terminal_pair(
    profile: AssayProfile,
    left_matches: Mapping[str, TerminalPrimerMatch],
    right_matches: Mapping[str, TerminalPrimerMatch],
) -> tuple[
    ProductDefinition,
    TerminalPrimerMatch,
    TerminalPrimerMatch,
    int,
] | None:
    best = None
    best_score = float("-inf")
    for product in profile.products:
        left = left_matches.get(product.forward_primer)
        right = right_matches.get(product.reverse_primer)
        if left is None or right is None or left.read_end0 > right.read_start0:
            continue
        observed_length = right.read_end0 - left.read_start0
        score = left.score + right.score
        if score > best_score:
            best_score = score
            best = (product, left, right, observed_length)
    return best


def _best_incompatible_pair(
    left_matches: Mapping[str, TerminalPrimerMatch],
    right_matches: Mapping[str, TerminalPrimerMatch],
    declared_pairs: set[tuple[str, str]],
) -> tuple[TerminalPrimerMatch, TerminalPrimerMatch] | None:
    best = None
    best_score = float("-inf")
    for left_id, left in left_matches.items():
        for right_id, right in right_matches.items():
            if left.read_end0 > right.read_start0:
                continue
            if (left_id, right_id) in declared_pairs:
                continue
            score = left.score + right.score
            if score > best_score:
                best_score = score
                best = (left, right)
    return best


def _noncomplete_observation(
    *,
    record: FastqRecord,
    sample: str,
    profile: AssayProfile,
    evaluated: Sequence[
        tuple[
            str,
            str,
            str,
            list[_EndpointCandidate],
            Mapping[str, TerminalPrimerMatch],
            Mapping[str, TerminalPrimerMatch],
        ]
    ],
    config: AdmissionConfig,
) -> MoleculeObservation:
    """Classify molecules lacking a valid declared product span."""

    unexpected_declared = []
    for orientation, _, _, _, left_matches, right_matches in evaluated:
        pair = _best_declared_terminal_pair(
            profile, left_matches, right_matches
        )
        if pair:
            product, left, right, observed_length = pair
            unexpected_declared.append(
                (
                    left.score + right.score,
                    orientation,
                    product,
                    left,
                    right,
                    observed_length,
                )
            )
    if unexpected_declared:
        _, orientation, product, left, right, observed_length = max(
            unexpected_declared, key=lambda item: item[0]
        )
        lower_tolerance = max(
            2, int(product.min_length * config.length_tolerance_fraction)
        )
        status = (
            "primer_dimer_candidate"
            if observed_length < product.min_length - lower_tolerance
            else "unexpected_length"
        )
        return MoleculeObservation(
            sample=sample,
            read_id=record.read_id,
            read_length_bp=len(record.sequence),
            status=status,
            orientation=orientation,
            left_primer_id=left.primer_id,
            right_primer_id=right.primer_id,
            left_overlap_bp=left.overlap_bp,
            right_overlap_bp=right.overlap_bp,
            left_mismatches=left.mismatches,
            right_mismatches=right.mismatches,
            trim_start0=left.read_start0,
            trim_end0=right.read_end0,
            observed_product_length_bp=observed_length,
            candidate_product_ids=(product.product_id,),
            reason=(
                f"declared primer pair span {observed_length} bp is outside "
                f"product range {product.min_length}-{product.max_length} bp"
            ),
        )

    declared_pairs = {
        (product.forward_primer, product.reverse_primer)
        for product in profile.products
    }
    incompatible = []
    for orientation, _, _, _, left_matches, right_matches in evaluated:
        pair = _best_incompatible_pair(
            left_matches, right_matches, declared_pairs
        )
        if pair:
            incompatible.append(
                (pair[0].score + pair[1].score, orientation, pair)
            )
    if incompatible:
        _, orientation, (left, right) = max(
            incompatible, key=lambda item: item[0]
        )
        return MoleculeObservation(
            sample=sample,
            read_id=record.read_id,
            read_length_bp=len(record.sequence),
            status="chimera_candidate",
            orientation=orientation,
            left_primer_id=left.primer_id,
            right_primer_id=right.primer_id,
            left_overlap_bp=left.overlap_bp,
            right_overlap_bp=right.overlap_bp,
            left_mismatches=left.mismatches,
            right_mismatches=right.mismatches,
            trim_start0=left.read_start0,
            trim_end0=right.read_end0,
            observed_product_length_bp=right.read_end0 - left.read_start0,
            reason="terminal primers form no declared assay product",
        )

    best_one = []
    for orientation, _, _, _, left_matches, right_matches in evaluated:
        for match in left_matches.values():
            best_one.append((match.score, orientation, "left", match))
        for match in right_matches.values():
            best_one.append((match.score, orientation, "right", match))
    if best_one:
        _, orientation, side, match = max(
            best_one, key=lambda item: item[0]
        )
        return MoleculeObservation(
            sample=sample,
            read_id=record.read_id,
            read_length_bp=len(record.sequence),
            status="one_ended",
            orientation=orientation,
            left_primer_id=match.primer_id if side == "left" else "",
            right_primer_id=match.primer_id if side == "right" else "",
            left_overlap_bp=match.overlap_bp if side == "left" else None,
            right_overlap_bp=match.overlap_bp if side == "right" else None,
            left_mismatches=match.mismatches if side == "left" else None,
            right_mismatches=match.mismatches if side == "right" else None,
            reason=f"only a {side}-terminal primer was recognised",
        )

    return MoleculeObservation(
        sample=sample,
        read_id=record.read_id,
        read_length_bp=len(record.sequence),
        status="off_target",
        orientation="unknown",
        reason="no declared terminal primer was recognised",
    )


def _complete_observation(
    *,
    record: FastqRecord,
    sample: str,
    profile: AssayProfile,
    compiled: CompiledAssay,
    config: AdmissionConfig,
    full: Sequence[tuple[float, str, str, str, _EndpointCandidate]],
    compiled_index: CompiledSequenceIndex | None,
) -> MoleculeObservation:
    _, orientation, oriented_sequence, _, top = full[0]

    product_scores: dict[str, float] = {}
    product_candidates: dict[str, _EndpointCandidate] = {}
    for score, candidate_orientation, _, _, candidate in full:
        if candidate_orientation != orientation:
            continue
        product_id = candidate.product.product_id
        if product_id not in product_scores or score > product_scores[product_id]:
            product_scores[product_id] = score
            product_candidates[product_id] = candidate

    candidate_ids = tuple(
        product_id
        for product_id, _ in sorted(
            product_scores.items(), key=lambda item: (-item[1], item[0])
        )
    )
    assigned_product = candidate_ids[0]
    if len(candidate_ids) > 1:
        endpoint_margin = (
            product_scores[candidate_ids[0]]
            - product_scores[candidate_ids[1]]
        )
        if endpoint_margin < 1.0:
            return MoleculeObservation(
                sample=sample,
                read_id=record.read_id,
                read_length_bp=len(record.sequence),
                status="complete",
                orientation=orientation,
                left_primer_id=top.left.primer_id,
                right_primer_id=top.right.primer_id,
                left_overlap_bp=top.left.overlap_bp,
                right_overlap_bp=top.right.overlap_bp,
                left_mismatches=top.left.mismatches,
                right_mismatches=top.right.mismatches,
                trim_start0=top.observed_start0,
                trim_end0=top.observed_end0,
                observed_product_length_bp=top.observed_length_bp,
                candidate_product_ids=candidate_ids,
                assignment_state="ambiguous_product_definition",
                reason=(
                    "multiple declared products have indistinguishable "
                    "endpoint evidence"
                ),
            )

    endpoint = product_candidates[assigned_product]
    observed = oriented_sequence[
        endpoint.observed_start0 : endpoint.observed_end0
    ]
    left_missing = endpoint.left.primer_start0
    reverse_primer_len = len(
        profile.primers[endpoint.product.reverse_primer].sequence
    )
    right_missing = reverse_primer_len - endpoint.right.primer_end0
    sequence_scores = score_compiled_sequences(
        observed,
        product_id=endpoint.product.product_id,
        left_missing_bp=left_missing,
        right_missing_bp=right_missing,
        products=compiled.products,
        config=config,
        compiled_index=compiled_index,
    )

    base_kwargs = dict(
        sample=sample,
        read_id=record.read_id,
        read_length_bp=len(record.sequence),
        status="complete",
        orientation=orientation,
        left_primer_id=endpoint.left.primer_id,
        right_primer_id=endpoint.right.primer_id,
        left_overlap_bp=endpoint.left.overlap_bp,
        right_overlap_bp=endpoint.right.overlap_bp,
        left_mismatches=endpoint.left.mismatches,
        right_mismatches=endpoint.right.mismatches,
        trim_start0=endpoint.observed_start0,
        trim_end0=endpoint.observed_end0,
        observed_product_length_bp=endpoint.observed_length_bp,
        candidate_product_ids=candidate_ids,
        assigned_product_id=assigned_product,
    )
    if not sequence_scores:
        return MoleculeObservation(
            **base_kwargs,
            assignment_state="no_compiled_sequence",
            reason="endpoint product has no compiled sequence candidate",
        )

    best = sequence_scores[0]
    second = sequence_scores[1] if len(sequence_scores) > 1 else None
    margin = second.edit_rate - best.edit_rate if second else None
    if best.edit_rate > config.max_sequence_edit_rate:
        state = "poor_sequence_fit"
        reason = (
            f"best compiled sequence edit rate {best.edit_rate:.4f} exceeds "
            f"{config.max_sequence_edit_rate:.4f}"
        )
        haplotypes = best.haplotype_ids
    elif (
        second
        and margin is not None
        and margin < config.sequence_assignment_margin
    ):
        state = "low_margin"
        near = [
            item
            for item in sequence_scores
            if item.edit_rate - best.edit_rate
            < config.sequence_assignment_margin
        ]
        haplotypes = tuple(
            sorted({hap for item in near for hap in item.haplotype_ids})
        )
        reason = "multiple compiled sequences have similar edit rates"
    elif len(best.haplotype_ids) > 1:
        state = "equivalent_haplotypes"
        haplotypes = best.haplotype_ids
        reason = (
            "multiple haplotypes compile to the same observed product sequence"
        )
    else:
        state = "unique_sequence"
        haplotypes = best.haplotype_ids
        reason = "one compiled product sequence is best supported"

    return MoleculeObservation(
        **base_kwargs,
        assignment_state=state,
        best_sequence_sha256=best.sequence_sha256,
        candidate_haplotypes=haplotypes,
        edit_distance=best.edit_distance,
        edit_rate=best.edit_rate,
        sequence_log_likelihood=best.log_likelihood,
        likelihood_margin=margin,
        reason=reason,
    )


PrimerCandidateProvider = Callable[[str, str], Sequence[str]]


def admit_molecule(
    record: FastqRecord,
    *,
    sample: str,
    profile: AssayProfile,
    compiled: CompiledAssay,
    config: AdmissionConfig | None = None,
    _primer_candidate_provider: PrimerCandidateProvider | None = None,
    _compiled_index: CompiledSequenceIndex | None = None,
) -> MoleculeObservation:
    """Classify one FASTQ molecule and score compatible compiled products."""

    config = config or AdmissionConfig()
    config.validate()
    if (
        compiled.assay_id != profile.assay_id
        or compiled.assay_version != profile.version
    ):
        raise MoleculeAdmissionError("compiled assay/profile mismatch")

    orientations = [
        ("forward", record.sequence, record.qualities),
        (
            "reverse",
            reverse_complement(record.sequence),
            _reverse_qualities(record.qualities),
        ),
    ]
    evaluated = []
    for orientation, sequence, qualities in orientations:
        if _primer_candidate_provider is None:
            left_ids = right_ids = None
        else:
            left_ids = _primer_candidate_provider(sequence, "left")
            right_ids = _primer_candidate_provider(sequence, "right")
        candidates, left_matches, right_matches = (
            _endpoint_candidates_for_orientation(
                profile,
                sequence,
                qualities,
                orientation,
                config,
                left_primer_ids=left_ids,
                right_primer_ids=right_ids,
            )
        )
        evaluated.append(
            (
                orientation,
                sequence,
                qualities,
                candidates,
                left_matches,
                right_matches,
            )
        )

    full = [
        (candidate.endpoint_score, orientation, sequence, qualities, candidate)
        for orientation, sequence, qualities, candidates, _, _ in evaluated
        for candidate in candidates
    ]
    full.sort(
        key=lambda item: (-item[0], item[4].product.product_id, item[1])
    )
    if not full:
        return _noncomplete_observation(
            record=record,
            sample=sample,
            profile=profile,
            evaluated=evaluated,
            config=config,
        )
    return _complete_observation(
        record=record,
        sample=sample,
        profile=profile,
        compiled=compiled,
        config=config,
        full=full,
        compiled_index=_compiled_index,
    )


class AdmissionEngine:
    """Reusable per-assay engine with primer-seed and product indexes."""

    def __init__(
        self,
        *,
        profile: AssayProfile,
        compiled: CompiledAssay,
        config: AdmissionConfig | None = None,
    ) -> None:
        self.profile = profile
        self.compiled = compiled
        self.config = config or AdmissionConfig()
        self.config.validate()
        if (
            compiled.assay_id != profile.assay_id
            or compiled.assay_version != profile.version
        ):
            raise MoleculeAdmissionError("compiled assay/profile mismatch")

        self._all_primers = tuple(sorted(profile.primers))
        self._max_primer_len = max(
            len(primer.sequence) for primer in profile.primers.values()
        )
        self._left_seed_index, self._left_unindexed = self._build_seed_index(
            {pid: primer.sequence for pid, primer in profile.primers.items()}
        )
        self._right_seed_index, self._right_unindexed = self._build_seed_index(
            {
                pid: reverse_complement(primer.sequence)
                for pid, primer in profile.primers.items()
            }
        )
        self._compiled_index = build_compiled_sequence_index(compiled.products)

    def _build_seed_index(
        self,
        patterns: Mapping[str, str],
    ) -> tuple[dict[str, set[str]], set[str]]:
        seed_length = self.config.primer_seed_length
        index: dict[str, set[str]] = defaultdict(set)
        unindexed: set[str] = set()
        for primer_id, pattern in patterns.items():
            concrete = normalise_dna(pattern)
            kmers = {
                concrete[start : start + seed_length]
                for start in range(0, len(concrete) - seed_length + 1)
                if set(concrete[start : start + seed_length]) <= set("ACGT")
            }
            if not kmers:
                unindexed.add(primer_id)
                continue
            for kmer in kmers:
                index[kmer].add(primer_id)
        return dict(index), unindexed

    def _candidate_primer_ids(
        self,
        sequence: str,
        side: str,
    ) -> tuple[str, ...]:
        search_width = min(
            len(sequence), self.config.end_search_bp + self._max_primer_len
        )
        window = (
            sequence[:search_width]
            if side == "left"
            else sequence[-search_width:]
        )
        index = (
            self._left_seed_index if side == "left" else self._right_seed_index
        )
        unindexed = (
            self._left_unindexed if side == "left" else self._right_unindexed
        )
        seed_length = self.config.primer_seed_length
        candidates: set[str] = set(unindexed)
        for start in range(0, len(window) - seed_length + 1):
            candidates.update(index.get(window[start : start + seed_length], ()))

        # Exact seeding is a performance shortcut, never a scientific gate.
        if candidates == unindexed:
            return self._all_primers
        return tuple(sorted(candidates))

    @staticmethod
    def _evidence_rank(value: MoleculeObservation) -> tuple[int, int, float]:
        status_rank = {
            "complete": 6,
            "primer_dimer_candidate": 5,
            "unexpected_length": 4,
            "chimera_candidate": 3,
            "one_ended": 2,
            "off_target": 1,
        }
        assignment_rank = {
            "unique_sequence": 5,
            "equivalent_haplotypes": 4,
            "low_margin": 3,
            "poor_sequence_fit": 2,
            "no_compiled_sequence": 1,
            "ambiguous_product_definition": 1,
            "": 0,
        }
        edit = value.edit_rate if value.edit_rate is not None else 1.0
        return (
            status_rank.get(value.status, 0),
            assignment_rank.get(value.assignment_state, 0),
            -edit,
        )

    def admit(self, record: FastqRecord, *, sample: str) -> MoleculeObservation:
        observation = admit_molecule(
            record,
            sample=sample,
            profile=self.profile,
            compiled=self.compiled,
            config=self.config,
            _primer_candidate_provider=self._candidate_primer_ids,
            _compiled_index=self._compiled_index,
        )

        # Retry unresolved reads exhaustively so a noisy true primer is not lost
        # merely because all of its short exact seeds contain an error.
        questionable = {
            "ambiguous_product_definition",
            "no_compiled_sequence",
            "poor_sequence_fit",
            "low_margin",
        }
        if (
            observation.status != "complete"
            or observation.assignment_state in questionable
        ):
            exhaustive = admit_molecule(
                record,
                sample=sample,
                profile=self.profile,
                compiled=self.compiled,
                config=self.config,
                _compiled_index=self._compiled_index,
            )
            if self._evidence_rank(exhaustive) > self._evidence_rank(observation):
                observation = exhaustive
        return observation


__all__ = [
    "AdmissionAccumulator",
    "AdmissionConfig",
    "AdmissionEngine",
    "CompiledAssay",
    "FastqRecord",
    "MoleculeAdmissionError",
    "MoleculeObservation",
    "OBSERVATION_FIELDS",
    "admit_molecule",
    "load_compiled_assay",
    "read_fastq",
    "summary_payload",
    "write_observation_tsv_header",
    "write_product_counts_tsv",
]
