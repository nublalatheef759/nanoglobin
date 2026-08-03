#!/usr/bin/env python3
"""Fit product efficiency, dropout and artifact terms from declared controls."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.calibration_fit import calibrate_model  # noqa: E402
from nanoglobin.calibration_io import (  # noqa: E402
    load_calibration_settings,
    product_rows,
    read_control_manifest,
    sample_product_rows,
    write_model_config,
    write_stratum_configs,
    write_tsv,
)
from nanoglobin.calibration_types import CalibrationError  # noqa: E402
from nanoglobin.genotype import (  # noqa: E402
    GenotypeModelError,
    load_compiled_space,
    load_config,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit relative PCR-product efficiencies, required-product dropout and "
            "unsupported-molecule mass from a truth-scoped control manifest."
        )
    )
    parser.add_argument("--compiled-json", required=True)
    parser.add_argument("--control-manifest", required=True)
    parser.add_argument("--base-model-config")
    parser.add_argument("--calibration-settings")
    parser.add_argument("--calibration-json", required=True)
    parser.add_argument("--calibrated-model-config", required=True)
    parser.add_argument("--product-parameters-tsv", required=True)
    parser.add_argument("--sample-product-fit-tsv", required=True)
    parser.add_argument("--stratum-config-dir")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        space = load_compiled_space(args.compiled_json)
        controls = read_control_manifest(args.control_manifest)
        base_config = load_config(args.base_model_config)
        settings = load_calibration_settings(args.calibration_settings)
        result = calibrate_model(
            space,
            controls,
            base_config=base_config,
            settings=settings,
        )
        calibration_path = Path(args.calibration_json)
        calibration_path.parent.mkdir(parents=True, exist_ok=True)
        calibration_path.write_text(
            json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        write_model_config(args.calibrated_model_config, result.calibrated_config)
        write_tsv(args.product_parameters_tsv, product_rows(result))
        write_tsv(args.sample_product_fit_tsv, sample_product_rows(result))
        stratum_paths = (
            write_stratum_configs(args.stratum_config_dir, result)
            if args.stratum_config_dir
            else []
        )
    except (
        CalibrationError,
        GenotypeModelError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        print(f"calibrate_amplicon_model: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "calibration_id": result.calibration_id,
                "calibration_status": result.calibration_status,
                "eligible_controls": sum(item.eligible for item in result.controls),
                "excluded_controls": sum(not item.eligible for item in result.controls),
                "product_count": sum(len(item.products) for item in result.components),
                "stratum_configs": [str(path) for path in stratum_paths],
                "warnings": list(result.warnings),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
