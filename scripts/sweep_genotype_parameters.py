#!/usr/bin/env python3
"""Evaluate genotype calls over a declared model-parameter grid."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.genotype import (  # noqa: E402
    GenotypeModelError,
    load_compiled_space,
    load_config,
    read_molecule_evidence,
)
from nanoglobin.stability import (  # noqa: E402
    StabilityError,
    load_stability_grid,
    run_stability_grid,
    stability_summary,
    write_stability_tsv,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the chromosome-haplotype posterior across a declared scalar "
            "parameter grid and report top-class/call-state stability."
        )
    )
    parser.add_argument("--compiled-json", required=True)
    parser.add_argument("--molecules-tsv", required=True)
    parser.add_argument("--model-config")
    parser.add_argument("--grid", required=True)
    parser.add_argument("--scenarios-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        space = load_compiled_space(args.compiled_json)
        evidence = read_molecule_evidence(args.molecules_tsv)
        config = load_config(args.model_config)
        grid = load_stability_grid(args.grid, config)
        scenarios = run_stability_grid(space, evidence, config, grid)
        summary = stability_summary(
            scenarios,
            base_config=config,
            grid=grid,
        )
        write_stability_tsv(args.scenarios_tsv, scenarios)
        output = Path(args.summary_json)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    except (
        StabilityError,
        GenotypeModelError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        print(f"sweep_genotype_parameters: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
