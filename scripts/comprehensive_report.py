#!/usr/bin/env python3
"""Build the per-variant report from an explicit sample manifest.

The previous implementation globbed ``results/`` and ``variants/``.  That made
workflow provenance depend on whatever happened to be present on disk and could
silently include simulation or stale samples in a clinical-facing report.  This
version accepts the sample IDs explicitly and derives only their declared files.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable, Iterator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.identify_sv import identify_sv  # noqa: E402


REPORT_COLUMNS = [
    "Sample",
    "Tool",
    "Chromosome",
    "Position",
    "Ref",
    "Alt",
    "HGVS",
    "Gene",
    "Consequence",
    "Genotype",
    "Quality",
]


def _require(path: Path, *, kind: str, sample: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"missing {kind} for sample {sample}: {path}")
    return path


def annotated_rows(sample: str, path: Path) -> Iterator[list[object]]:
    with _require(path, kind="annotated small-variant CSV", sample=sample).open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        reader = csv.DictReader(handle)
        required = {
            "Sample",
            "Chromosome",
            "Position",
            "Ref",
            "Alt",
            "HGVS",
            "Gene",
            "Consequence",
            "Genotype",
            "Quality",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            if row["Sample"] and row["Sample"] != sample:
                raise ValueError(
                    f"{path} contains sample {row['Sample']!r}, expected {sample!r}"
                )
            if row["HGVS"] == "non-globin":
                continue
            yield [
                sample,
                "Clair3",
                row["Chromosome"],
                row["Position"],
                row["Ref"],
                row["Alt"],
                row["HGVS"],
                row["Gene"],
                row["Consequence"],
                row["Genotype"],
                row["Quality"],
            ]


def _parse_info(info_text: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for field in info_text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            parsed[key] = value
    return parsed


def _in_globin_region(chrom: str, pos: int) -> tuple[bool, str]:
    if chrom == "chr11" and 5_225_000 <= pos <= 5_310_000:
        return True, "HBB"
    if chrom == "chr16" and 165_000 <= pos <= 180_000:
        return True, "HBA"
    return False, ""


def sv_rows(sample: str, path: Path, tool: str) -> Iterator[list[object]]:
    with _require(path, kind=f"{tool} VCF", sample=sample).open(
        encoding="utf-8"
    ) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                raise ValueError(f"malformed VCF record in {path}:{line_number}")
            chrom = parts[0]
            pos = int(parts[1])
            in_region, gene = _in_globin_region(chrom, pos)
            if not in_region:
                continue
            info = _parse_info(parts[7])
            svtype = info.get("SVTYPE", "")
            svlen = info.get("SVLEN", "")
            if svtype in {"INV", "BND"} or not svtype:
                continue
            genotype = parts[9].split(":", 1)[0] if len(parts) > 9 else ""
            name = identify_sv(chrom, str(pos), svtype, svlen)
            ref_display = parts[3][:20] + "..." if len(parts[3]) > 20 else parts[3]
            yield [
                sample,
                tool,
                chrom,
                pos,
                ref_display,
                svtype,
                name,
                gene,
                f"SV_{svtype}_{svlen}bp",
                genotype,
                parts[5],
            ]


def coverage_rows(sample: str, path: Path) -> Iterator[list[object]]:
    with _require(path, kind="coverage summary", sample=sample).open(
        encoding="utf-8"
    ) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                raise ValueError(f"malformed coverage row in {path}:{line_number}")
            yield [sample, "Coverage", "", "", "", "", "", parts[0], parts[1], "", ""]


def build_report(
    output: Path,
    samples: Iterable[str],
    *,
    results_dir: Path = Path("results"),
    variants_dir: Path = Path("variants"),
) -> None:
    sample_list = list(samples)
    if not sample_list:
        raise ValueError("at least one sample is required")
    if len(set(sample_list)) != len(sample_list):
        raise ValueError("sample list contains duplicates")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(REPORT_COLUMNS)
        for sample in sample_list:
            writer.writerows(annotated_rows(sample, results_dir / f"{sample}.annotated.csv"))
            sample_dir = variants_dir / sample
            writer.writerows(sv_rows(sample, sample_dir / f"{sample}.sniffles.vcf", "Sniffles"))
            writer.writerows(sv_rows(sample, sample_dir / f"{sample}.cutesv.vcf", "CuteSV"))
            writer.writerows(coverage_rows(sample, sample_dir / f"{sample}.coverage.tsv"))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a comprehensive report for exactly the declared samples."
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample", action="append", dest="samples", required=True)
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--variants-dir", default="variants")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        build_report(
            Path(args.output),
            args.samples,
            results_dir=Path(args.results_dir),
            variants_dir=Path(args.variants_dir),
        )
    except (OSError, ValueError) as exc:
        print(f"comprehensive_report: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
