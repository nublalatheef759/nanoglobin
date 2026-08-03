#!/usr/bin/env python3
"""Infer a diploid chromosome-haplotype posterior from molecule evidence."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.calibration_io import load_calibration_provenance  # noqa: E402
from nanoglobin.calibration_types import CalibrationError  # noqa: E402
from nanoglobin.genotype import (  # noqa: E402
    GENOTYPE_KERNEL,
    GenotypeModelError,
    analysis_payload,
    load_compiled_space,
    load_config,
    make_call,
    read_molecule_evidence,
    score_genotypes,
    write_scores_tsv,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score diploid chromosome-haplotype pairs from compiled assay "
            "products and primer-aware molecule evidence."
        )
    )
    parser.add_argument("--compiled-json", required=True)
    parser.add_argument("--molecules-tsv", required=True)
    parser.add_argument("--model-config")
    parser.add_argument(
        "--calibration-json",
        help=(
            "Optional provenance emitted by calibrate_amplicon_model.py. The "
            "embedded calibrated config must exactly match --model-config."
        ),
    )
    parser.add_argument("--posteriors-tsv", required=True)
    parser.add_argument("--call-json", required=True)
    parser.add_argument("--summary-json", required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        space = load_compiled_space(args.compiled_json)
        evidence = read_molecule_evidence(args.molecules_tsv)
        config = load_config(args.model_config)
        calibration = load_calibration_provenance(
            args.calibration_json,
            space,
            config,
        )
        scores, useful_effective_reads = score_genotypes(space, evidence, config)
        call = make_call(scores, evidence, useful_effective_reads, config)
        payload = analysis_payload(
            space=space,
            config=config,
            scores=scores,
            call=call,
        )
        if calibration is not None:
            payload["calibration"] = calibration
            payload["model"]["calibration_status"] = calibration[
                "calibration_status"
            ]
        write_scores_tsv(args.posteriors_tsv, scores)
        call_path = Path(args.call_json)
        call_path.parent.mkdir(parents=True, exist_ok=True)
        call_path.write_text(
            json.dumps(call.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        summary_path = Path(args.summary_json)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    except (
        CalibrationError,
        GenotypeModelError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        print(f"genotype_amplicons: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                **call.to_dict(),
                "kernel": GENOTYPE_KERNEL,
                "calibration_id": (
                    calibration["calibration_id"] if calibration else None
                ),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
