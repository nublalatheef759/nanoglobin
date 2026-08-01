"""Quality-aware terminal primer recognition for ONT amplicon reads."""

from __future__ import annotations

from typing import Iterator
import math

from .assay import IUPAC_BASES, normalise_dna
from .molecule_types import (
    AdmissionConfig,
    MoleculeAdmissionError,
    TerminalPrimerMatch,
)


def _base_compatible(observed: str, pattern: str) -> bool:
    # An observed N is missing evidence, not a match to every primer base.
    if observed not in "ACGT":
        return False
    return bool(IUPAC_BASES[observed] & IUPAC_BASES[pattern])


def _quality_log_likelihood(observed: str, pattern: str, qualities: str) -> float:
    value = 0.0
    for base, expected, char in zip(observed, pattern, qualities):
        q = max(0, min(60, ord(char) - 33))
        error = min(0.25, max(1e-6, 10 ** (-q / 10)))
        if _base_compatible(base, expected):
            value += math.log(max(1e-12, 1 - error))
        elif base == "N":
            value += math.log(0.25)
        else:
            value += math.log(max(1e-12, error / 3))
    return value


def _terminal_candidates(
    sequence: str,
    primer: str,
    *,
    side: str,
    config: AdmissionConfig,
) -> Iterator[tuple[int, int, int, int, int]]:
    """Yield read/primer spans plus terminal offset for one read side."""

    primer_len = len(primer)
    min_overlap = min(config.min_primer_overlap, primer_len)
    read_len = len(sequence)
    if read_len < min_overlap:
        return

    seen: set[tuple[int, int, int, int]] = set()
    if side == "left":
        max_start = min(config.end_search_bp, read_len - min_overlap)
        for read_start in range(max_start + 1):
            overlap = min(primer_len, read_len - read_start)
            if overlap < min_overlap:
                continue
            key = (read_start, read_start + overlap, 0, overlap)
            if key not in seen:
                seen.add(key)
                yield (*key, read_start)
        # The read starts after the primer began: its prefix is a primer suffix.
        for primer_start in range(1, primer_len - min_overlap + 1):
            overlap = min(primer_len - primer_start, read_len)
            if overlap < min_overlap:
                continue
            key = (0, overlap, primer_start, primer_start + overlap)
            if key not in seen:
                seen.add(key)
                yield (*key, 0)
    elif side == "right":
        max_offset = min(config.end_search_bp, read_len - min_overlap)
        for offset in range(max_offset + 1):
            read_end = read_len - offset
            if read_end < primer_len:
                continue
            key = (read_end - primer_len, read_end, 0, primer_len)
            if key not in seen:
                seen.add(key)
                yield (*key, offset)
        # The read ends before the primer is complete: suffix is primer prefix.
        for primer_end in range(min_overlap, primer_len):
            if read_len < primer_end:
                continue
            key = (read_len - primer_end, read_len, 0, primer_end)
            if key not in seen:
                seen.add(key)
                yield (*key, 0)
    else:
        raise ValueError("side must be left or right")


def best_terminal_primer_match(
    sequence: str,
    qualities: str,
    *,
    primer_id: str,
    primer_sequence: str,
    side: str,
    config: AdmissionConfig,
) -> TerminalPrimerMatch | None:
    """Return the best accepted quality-aware terminal primer match."""

    config.validate()
    sequence_n = normalise_dna(sequence)
    primer_n = normalise_dna(primer_sequence, label=f"primer {primer_id}")
    if len(sequence_n) != len(qualities):
        raise MoleculeAdmissionError("sequence and quality lengths differ")

    best: TerminalPrimerMatch | None = None
    for read_start, read_end, primer_start, primer_end, offset in _terminal_candidates(
        sequence_n, primer_n, side=side, config=config
    ):
        observed = sequence_n[read_start:read_end]
        pattern = primer_n[primer_start:primer_end]
        overlap = len(observed)
        mismatches = sum(
            0 if _base_compatible(base, expected) else 1
            for base, expected in zip(observed, pattern)
        )
        allowed_by_rate = int(math.floor(overlap * config.max_primer_error_rate))
        allowed = min(config.max_primer_mismatches, allowed_by_rate)
        if mismatches > allowed:
            continue
        quality_ll = _quality_log_likelihood(
            observed, pattern, qualities[read_start:read_end]
        )
        clipped = primer_start + (len(primer_n) - primer_end)
        score = (
            overlap
            - 4.0 * mismatches
            - config.clipped_primer_penalty * clipped
            - config.terminal_offset_penalty * offset
            + 0.05 * quality_ll
        )
        candidate = TerminalPrimerMatch(
            primer_id=primer_id,
            side=side,
            read_start0=read_start,
            read_end0=read_end,
            primer_start0=primer_start,
            primer_end0=primer_end,
            primer_length_bp=len(primer_n),
            overlap_bp=overlap,
            mismatches=mismatches,
            offset_bp=offset,
            quality_log_likelihood=quality_ll,
            score=score,
        )
        if best is None or (
            candidate.score,
            candidate.overlap_bp,
            -candidate.mismatches,
            -candidate.offset_bp,
        ) > (
            best.score,
            best.overlap_bp,
            -best.mismatches,
            -best.offset_bp,
        ):
            best = candidate
    return best
