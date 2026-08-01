#!/usr/bin/env python3
"""Build an HBA family-coordinate and copy-marker map from candidate sequences."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.family_map import (  # noqa: E402
    FamilyMapError,
    build_coordinate_map,
    marker_rows,
    read_fasta,
    summary_payload,
    write_map_json,
    write_markers_tsv,
    write_positions_tsv,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Align HBA1/HBA2/hybrid/candidate-copy sequences to one family "
            "coordinate anchor and report copy-informative markers."
        )
    )
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--anchor-id", required=True)
    parser.add_argument("--positions-tsv", required=True)
    parser.add_argument("--markers-tsv", required=True)
    parser.add_argument("--map-json", required=True)
    parser.add_argument("--summary-json", required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        records = read_fasta(args.sequences)
        anchor, positions = build_coordinate_map(records, anchor_id=args.anchor_id)
        markers = marker_rows(records, positions, anchor_id=args.anchor_id)
        summary = summary_payload(
            anchor=anchor,
            records=records,
            observations=positions,
            markers=markers,
        )
        write_positions_tsv(args.positions_tsv, positions)
        write_markers_tsv(args.markers_tsv, markers)
        write_map_json(
            args.map_json,
            anchor=anchor,
            records=records,
            observations=positions,
            markers=markers,
        )
        summary_path = Path(args.summary_json)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    except (FamilyMapError, OSError, ValueError) as exc:
        print(f"build_hba_family_map: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
