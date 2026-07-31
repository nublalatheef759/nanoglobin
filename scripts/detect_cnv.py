#!/usr/bin/env python3
"""Detect local coverage departures without over-naming HBA copy gains.

Coverage is an evidence producer.  A gain-shaped interval can be compared with
IthaCNVs as a *candidate class*, but reciprocal overlap cannot establish copy
order, junction sequence, or a clinical allele such as alpha-alpha-alpha
anti-3.7.  The output therefore keeps candidate naming separate from the call.

The command accepts explicit ``SAMPLE=PATH`` inputs.  It never scans the
filesystem, which prevents stale or simulated samples from entering reports.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.identify_sv import identify_sv  # noqa: E402


MIN_RUN = 3
GAIN_FACTOR = 1.5
LOSS_FACTOR = 0.5
BASELINE = 1.0


def load_bins(path: str | Path, region: str = "HBA") -> list[tuple[int, int, float]]:
    rows: list[tuple[int, int, float]] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"region", "bin_start", "bin_end", "ratio_to_HBB"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            if row["region"] != region:
                continue
            rows.append((int(row["bin_start"]), int(row["bin_end"]), float(row["ratio_to_HBB"])))
    rows.sort()
    return rows


def find_runs(
    rows: Sequence[tuple[int, int, float]], baseline: float = BASELINE
) -> list[tuple[str, int, int, float, int | None]]:
    calls: list[tuple[str, int, int, float, int | None]] = []

    def classify(ratio: float) -> str | None:
        if baseline <= 0 or ratio != ratio:
            return None
        fold = ratio / baseline
        if fold >= GAIN_FACTOR:
            return "gain"
        if fold <= LOSS_FACTOR:
            return "loss"
        return None

    i = 0
    while i < len(rows):
        kind = classify(rows[i][2])
        if kind is None:
            i += 1
            continue
        j = i
        while j < len(rows) and classify(rows[j][2]) == kind:
            j += 1
        if j - i >= MIN_RUN:
            segment = rows[i:j]
            start = segment[0][0]
            end = segment[-1][1]
            median_ratio = statistics.median(ratio for _, _, ratio in segment)
            copies = round(2 * median_ratio / baseline) if baseline else None
            calls.append((kind, start, end, median_ratio / baseline, copies))
        i = j
    return calls


def call_file(path: str | Path) -> tuple[list[tuple[str, int, int, float, int | None]], float] | None:
    rows = load_bins(path)
    if not rows:
        return None
    return find_runs(rows, BASELINE), BASELINE


def _parse_bin_spec(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--bin must be SAMPLE=PATH")
    sample, path = value.split("=", 1)
    sample = sample.strip()
    if not sample or not path.strip():
        raise argparse.ArgumentTypeError("--bin must contain a non-empty sample and path")
    return sample, Path(path)


def write_cnv_calls(
    out_path: str | Path,
    bin_specs: Iterable[tuple[str, Path]],
) -> None:
    specs = list(bin_specs)
    if not specs:
        raise ValueError("at least one explicit coverage-bin input is required")
    sample_names = [sample for sample, _ in specs]
    if len(set(sample_names)) != len(sample_names):
        raise ValueError("duplicate sample in --bin inputs")

    output = Path(out_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Sample",
                "Type",
                "Region",
                "Fold",
                "Estimated_copies",
                "Name",
                "Naming_status",
                "Confidence",
            ]
        )
        for sample, path in specs:
            if not path.is_file():
                raise FileNotFoundError(f"missing coverage bins for {sample}: {path}")
            result = call_file(path)
            if result is None:
                continue
            calls, _baseline = result
            for kind, start, end, fold, copies in calls:
                region = f"chr16:{start}-{end}"
                name = ""
                naming_status = "not_applicable"
                if kind == "gain":
                    hit = identify_sv("chr16", start, "DUP", end - start + 1)
                    if hit and not hit.startswith("unknown"):
                        name = hit
                        naming_status = "catalogue_candidate_only"
                    else:
                        naming_status = "uncatalogued_gain_candidate"
                elif kind == "loss":
                    naming_status = "coverage_loss_not_typed"
                confidence = "moderate" if kind == "gain" else "screening"
                writer.writerow(
                    [
                        sample,
                        kind,
                        region,
                        f"{fold:.2f}",
                        "" if copies is None else copies,
                        name,
                        naming_status,
                        confidence,
                    ]
                )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--bin",
        action="append",
        required=True,
        type=_parse_bin_spec,
        metavar="SAMPLE=PATH",
        help="explicit per-sample coverage-bin TSV; repeat for each sample",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        write_cnv_calls(args.output, args.bin)
    except (OSError, ValueError) as exc:
        print(f"detect_cnv: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
