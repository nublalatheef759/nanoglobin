"""Strict streaming FASTQ input for NanoGlobin molecule evidence."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator
import gzip

from .assay import AssayProfileError, normalise_dna
from .molecule_types import FastqRecord, MoleculeAdmissionError


def _open_text(path: str | Path):
    value = str(path)
    if value.endswith(".gz"):
        return gzip.open(value, "rt", encoding="utf-8")
    return open(value, "r", encoding="utf-8")


def read_fastq(path: str | Path) -> Iterator[FastqRecord]:
    """Stream four-line FASTQ or FASTQ.GZ records with strict validation."""

    with _open_text(path) as handle:
        line_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            line_number += 1
            sequence = handle.readline()
            plus = handle.readline()
            qualities = handle.readline()
            if not sequence or not plus or not qualities:
                raise MoleculeAdmissionError(
                    f"truncated FASTQ record beginning at line {line_number} in {path}"
                )
            line_number += 3
            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            qualities = qualities.rstrip("\r\n")
            if not header.startswith("@"):
                raise MoleculeAdmissionError(
                    f"FASTQ header at line {line_number - 3} does not begin with @"
                )
            if not plus.startswith("+"):
                raise MoleculeAdmissionError(
                    f"FASTQ separator at line {line_number - 1} does not begin with +"
                )
            read_id = header[1:].split(maxsplit=1)[0]
            if not read_id:
                raise MoleculeAdmissionError(
                    f"empty FASTQ read identifier at line {line_number - 3}"
                )
            try:
                sequence_n = normalise_dna(sequence, label=f"FASTQ read {read_id}")
            except AssayProfileError as exc:
                raise MoleculeAdmissionError(str(exc)) from exc
            if len(sequence_n) != len(qualities):
                raise MoleculeAdmissionError(
                    f"FASTQ read {read_id} has sequence length {len(sequence_n)} "
                    f"but quality length {len(qualities)}"
                )
            if any(ord(char) < 33 or ord(char) > 126 for char in qualities):
                raise MoleculeAdmissionError(
                    f"FASTQ read {read_id} has invalid quality text"
                )
            yield FastqRecord(
                read_id=read_id, sequence=sequence_n, qualities=qualities
            )
