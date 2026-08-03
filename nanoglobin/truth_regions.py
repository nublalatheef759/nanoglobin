"""Audit whether declared globin targets are inside benchmark confidence regions."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping
import csv


class TruthRegionError(ValueError):
    """Raised when a BED or requested target region is invalid."""


@dataclass(frozen=True)
class TargetRegion:
    name: str
    chrom: str
    start0: int
    end0: int

    @property
    def length_bp(self) -> int:
        return self.end0 - self.start0


@dataclass(frozen=True)
class RegionAudit:
    name: str
    chrom: str
    start1: int
    end1: int
    length_bp: int
    covered_bp: int
    uncovered_bp: int
    coverage_fraction: float
    fully_covered: bool
    overlapping_segments: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def parse_target(value: str) -> TargetRegion:
    """Parse ``NAME=chrom:start-end`` with 1-based inclusive coordinates."""

    try:
        name, interval = value.split("=", 1)
        chrom, coordinates = interval.rsplit(":", 1)
        start_text, end_text = coordinates.replace(",", "").split("-", 1)
        start1 = int(start_text)
        end1 = int(end_text)
    except (ValueError, TypeError) as exc:
        raise TruthRegionError(
            f"invalid target {value!r}; expected NAME=chrom:start-end"
        ) from exc
    name = name.strip()
    chrom = chrom.strip()
    if not name or not chrom or start1 < 1 or end1 < start1:
        raise TruthRegionError(
            f"invalid target {value!r}; require non-empty name/chrom and 1 <= start <= end"
        )
    return TargetRegion(name=name, chrom=chrom, start0=start1 - 1, end0=end1)


def read_bed(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """Read and merge a standard 0-based half-open BED confidence file."""

    intervals: dict[str, list[tuple[int, int]]] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                raise TruthRegionError(
                    f"BED line {line_number} has fewer than three fields"
                )
            chrom = fields[0]
            try:
                start0 = int(fields[1])
                end0 = int(fields[2])
            except ValueError as exc:
                raise TruthRegionError(
                    f"BED line {line_number} has non-integer coordinates"
                ) from exc
            if start0 < 0 or end0 <= start0:
                raise TruthRegionError(
                    f"BED line {line_number} has invalid interval {start0}-{end0}"
                )
            intervals.setdefault(chrom, []).append((start0, end0))

    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, values in intervals.items():
        values.sort()
        chrom_merged: list[list[int]] = []
        for start0, end0 in values:
            if not chrom_merged or start0 > chrom_merged[-1][1]:
                chrom_merged.append([start0, end0])
            else:
                chrom_merged[-1][1] = max(chrom_merged[-1][1], end0)
        merged[chrom] = [(start0, end0) for start0, end0 in chrom_merged]
    return merged


def audit_regions(
    confidence: Mapping[str, Iterable[tuple[int, int]]],
    targets: Iterable[TargetRegion],
) -> list[RegionAudit]:
    output: list[RegionAudit] = []
    seen_names: set[str] = set()
    for target in targets:
        if target.name in seen_names:
            raise TruthRegionError(f"duplicate target name {target.name!r}")
        seen_names.add(target.name)
        covered = 0
        segments = 0
        for start0, end0 in confidence.get(target.chrom, []):
            if end0 <= target.start0:
                continue
            if start0 >= target.end0:
                break
            overlap = min(end0, target.end0) - max(start0, target.start0)
            if overlap > 0:
                covered += overlap
                segments += 1
        length = target.length_bp
        output.append(
            RegionAudit(
                name=target.name,
                chrom=target.chrom,
                start1=target.start0 + 1,
                end1=target.end0,
                length_bp=length,
                covered_bp=covered,
                uncovered_bp=length - covered,
                coverage_fraction=covered / length,
                fully_covered=covered == length,
                overlapping_segments=segments,
            )
        )
    return output


def write_audit_tsv(path: str | Path, audits: Iterable[RegionAudit]) -> None:
    rows = [item.to_dict() for item in audits]
    if not rows:
        raise TruthRegionError("at least one target is required")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
