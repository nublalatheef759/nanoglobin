#!/usr/bin/env python3
"""Benchmark and cross-check Python and Cython genotype scoring backends."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import mean, median, pstdev
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanoglobin.fast_genotype import (  # noqa: E402
    cython_backend_available,
    score_genotypes_auto,
)
from nanoglobin.genotype import (  # noqa: E402
    GenotypeModelError,
    load_compiled_space,
    load_config,
    read_molecule_evidence,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify numerical parity and measure Python/Cython chromosome-haplotype "
            "posterior runtimes on one real or synthetic molecule table."
        )
    )
    parser.add_argument("--compiled-json", required=True)
    parser.add_argument("--molecules-tsv", required=True)
    parser.add_argument("--model-config")
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-10)
    parser.add_argument("--relative-tolerance", type=float, default=1e-10)
    parser.add_argument("--output-json")
    return parser.parse_args(argv)


def _close(left: float, right: float, *, absolute: float, relative: float) -> bool:
    difference = abs(left - right)
    return difference <= absolute + relative * max(abs(left), abs(right))


def _assert_parity(reference, candidate, *, absolute: float, relative: float) -> None:
    if [item.class_id for item in reference] != [item.class_id for item in candidate]:
        raise GenotypeModelError("backend class ranking differs")
    for left, right in zip(reference, candidate):
        for field in (
            "log_read_likelihood",
            "log_count_likelihood",
            "log_dropout_likelihood",
            "log_prior",
            "log_score",
            "posterior",
        ):
            left_value = float(getattr(left, field))
            right_value = float(getattr(right, field))
            if not _close(
                left_value,
                right_value,
                absolute=absolute,
                relative=relative,
            ):
                raise GenotypeModelError(
                    f"backend mismatch for {left.class_id} {field}: "
                    f"python={left_value:.17g}, cython={right_value:.17g}"
                )


def _time_backend(space, evidence, config, backend: str, warmups: int, repeats: int):
    for _ in range(warmups):
        score_genotypes_auto(space, evidence, config, backend=backend)
    durations = []
    backend_used = ""
    for _ in range(repeats):
        start = perf_counter()
        _scores, _effective, backend_used = score_genotypes_auto(
            space, evidence, config, backend=backend
        )
        durations.append(perf_counter() - start)
    return {
        "backend_requested": backend,
        "backend_used": backend_used,
        "repeats": repeats,
        "seconds_min": min(durations),
        "seconds_median": median(durations),
        "seconds_mean": mean(durations),
        "seconds_max": max(durations),
        "seconds_population_sd": pstdev(durations),
        "all_seconds": durations,
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.warmups < 0 or args.repeats < 1:
        print("benchmark_genotype_backend: warmups >= 0 and repeats >= 1 required", file=sys.stderr)
        return 2
    try:
        space = load_compiled_space(args.compiled_json)
        evidence = read_molecule_evidence(args.molecules_tsv)
        config = load_config(args.model_config)
        reference, _effective, _ = score_genotypes_auto(
            space, evidence, config, backend="python"
        )
        results = [
            _time_backend(
                space,
                evidence,
                config,
                "python",
                args.warmups,
                args.repeats,
            )
        ]
        parity = "not_run"
        if cython_backend_available() and len(space.haplotype_ids) <= 64:
            candidate, _effective, _ = score_genotypes_auto(
                space, evidence, config, backend="cython"
            )
            _assert_parity(
                reference,
                candidate,
                absolute=args.absolute_tolerance,
                relative=args.relative_tolerance,
            )
            parity = "passed"
            results.append(
                _time_backend(
                    space,
                    evidence,
                    config,
                    "cython",
                    args.warmups,
                    args.repeats,
                )
            )
    except (GenotypeModelError, OSError, ValueError) as exc:
        print(f"benchmark_genotype_backend: {exc}", file=sys.stderr)
        return 2

    by_backend = {item["backend_requested"]: item for item in results}
    speedup = None
    if "cython" in by_backend:
        speedup = (
            by_backend["python"]["seconds_median"]
            / by_backend["cython"]["seconds_median"]
        )
    payload = {
        "schema_version": 1,
        "candidate_haplotype_count": len(space.haplotype_ids),
        "observable_class_count": len(reference),
        "molecule_count": len(evidence),
        "cython_extension_available": cython_backend_available(),
        "numerical_parity": parity,
        "median_speedup": speedup,
        "backends": results,
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output_json:
        output = Path(args.output_json)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
