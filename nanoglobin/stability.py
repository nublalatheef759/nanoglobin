"""Parameter-grid stability analysis for the chromosome-haplotype posterior."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, replace
from hashlib import sha256
from itertools import product
from pathlib import Path
from statistics import median
from typing import Any, Mapping, Sequence
import csv
import json

from .fast_genotype import BackendName, score_genotypes_auto
from .genotype import (
    CompiledGenotypeSpace,
    GenotypeConfig,
    GenotypeModelError,
    MoleculeEvidence,
    make_call,
)


class StabilityError(ValueError):
    """Raised when a sensitivity grid or result is invalid."""


_MAPPING_FIELDS = {"assignment_weights", "product_efficiencies", "haplotype_priors"}


@dataclass(frozen=True)
class StabilityGrid:
    axes: Mapping[str, tuple[int | float, ...]]
    max_scenarios: int = 1_000

    def validate(self, base_config: GenotypeConfig) -> None:
        if self.max_scenarios < 1:
            raise StabilityError("max_scenarios must be >= 1")
        if not self.axes:
            raise StabilityError("stability grid requires at least one axis")
        known = set(GenotypeConfig.__dataclass_fields__)
        unknown = sorted(set(self.axes) - known)
        if unknown:
            raise StabilityError(f"unknown genotype axis/axes: {', '.join(unknown)}")
        forbidden = sorted(set(self.axes) & _MAPPING_FIELDS)
        if forbidden:
            raise StabilityError(
                "mapping-valued axes require an explicit profile, not scalar grid "
                f"expansion: {', '.join(forbidden)}"
            )
        scenario_count = 1
        for name, values in self.axes.items():
            if not values:
                raise StabilityError(f"axis {name!r} has no values")
            current = getattr(base_config, name)
            for value in values:
                if isinstance(current, int) and not isinstance(current, bool):
                    if isinstance(value, bool) or int(value) != value:
                        raise StabilityError(f"axis {name!r} requires integer values")
                elif not isinstance(value, (int, float)) or isinstance(value, bool):
                    raise StabilityError(f"axis {name!r} requires numeric values")
            scenario_count *= len(values)
        if scenario_count > self.max_scenarios:
            raise StabilityError(
                f"grid expands to {scenario_count} scenarios; maximum is "
                f"{self.max_scenarios}"
            )
        # Validate every one-dimensional setting before the full Cartesian run so
        # errors identify the offending axis/value directly.
        for name, values in self.axes.items():
            current = getattr(base_config, name)
            for value in values:
                coerced: int | float
                if isinstance(current, int) and not isinstance(current, bool):
                    coerced = int(value)
                else:
                    coerced = float(value)
                try:
                    replace(base_config, **{name: coerced}).validate()
                except (GenotypeModelError, ValueError) as exc:
                    raise StabilityError(
                        f"invalid value {value!r} for axis {name!r}: {exc}"
                    ) from exc


@dataclass(frozen=True)
class StabilityScenario:
    scenario_id: str
    parameters: Mapping[str, int | float]
    call_state: str
    call_reason: str
    top_class_id: str
    top_posterior: float
    second_posterior: float | None
    posterior_odds: float | None
    genotype_pairs: tuple[tuple[str, str], ...]
    assigned_reads: int
    unresolved_fraction: float
    execution_backend: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "scenario_id": self.scenario_id,
            **dict(self.parameters),
            "call_state": self.call_state,
            "call_reason": self.call_reason,
            "top_class_id": self.top_class_id,
            "top_posterior": self.top_posterior,
            "second_posterior": self.second_posterior,
            "posterior_odds": self.posterior_odds,
            "genotype_pairs": ";".join(
                f"{left}/{right}" for left, right in self.genotype_pairs
            ),
            "assigned_reads": self.assigned_reads,
            "unresolved_fraction": self.unresolved_fraction,
            "execution_backend": self.execution_backend,
        }


def load_stability_grid(
    path: str | Path, base_config: GenotypeConfig
) -> StabilityGrid:
    profile = Path(path)
    text = profile.read_text(encoding="utf-8")
    if profile.suffix.lower() == ".json":
        raw = json.loads(text)
    else:
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("PyYAML is required to load a stability grid") from exc
        raw = yaml.safe_load(text)
    if not isinstance(raw, Mapping):
        raise StabilityError("stability profile must be a mapping")
    axes_raw = raw.get("axes")
    if not isinstance(axes_raw, Mapping):
        raise StabilityError("stability profile requires an axes mapping")
    axes: dict[str, tuple[int | float, ...]] = {}
    for name, values in axes_raw.items():
        if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
            raise StabilityError(f"axis {name!r} must be a list")
        axes[str(name)] = tuple(values)
    grid = StabilityGrid(
        axes=axes,
        max_scenarios=int(raw.get("max_scenarios", 1_000)),
    )
    grid.validate(base_config)
    return grid


def _scenarios(
    grid: StabilityGrid, base_config: GenotypeConfig
) -> list[tuple[dict[str, int | float], GenotypeConfig]]:
    names = sorted(grid.axes)
    output: list[tuple[dict[str, int | float], GenotypeConfig]] = []
    for values in product(*(grid.axes[name] for name in names)):
        parameters: dict[str, int | float] = {}
        updates: dict[str, int | float] = {}
        for name, value in zip(names, values):
            current = getattr(base_config, name)
            coerced: int | float
            if isinstance(current, int) and not isinstance(current, bool):
                coerced = int(value)
            else:
                coerced = float(value)
            parameters[name] = coerced
            updates[name] = coerced
        config = replace(base_config, **updates)
        config.validate()
        output.append((parameters, config))
    return output


def run_stability_grid(
    space: CompiledGenotypeSpace,
    evidence: Sequence[MoleculeEvidence],
    base_config: GenotypeConfig,
    grid: StabilityGrid,
    *,
    backend: BackendName = "auto",
) -> list[StabilityScenario]:
    grid.validate(base_config)
    output: list[StabilityScenario] = []
    for parameters, config in _scenarios(grid, base_config):
        canonical = json.dumps(parameters, sort_keys=True, separators=(",", ":"))
        scenario_id = "SG_" + sha256(canonical.encode("utf-8")).hexdigest()[:12]
        scores, effective, backend_used = score_genotypes_auto(
            space,
            evidence,
            config,
            backend=backend,
        )
        call = make_call(scores, evidence, effective, config)
        if call.top_class is None:
            raise StabilityError(f"scenario {scenario_id} produced no top class")
        second_posterior = (
            call.second_class.posterior if call.second_class is not None else None
        )
        odds = (
            call.top_class.posterior / max(second_posterior, 1e-300)
            if second_posterior is not None
            else None
        )
        output.append(
            StabilityScenario(
                scenario_id=scenario_id,
                parameters=parameters,
                call_state=call.call_state,
                call_reason=call.reason,
                top_class_id=call.top_class.class_id,
                top_posterior=call.top_class.posterior,
                second_posterior=second_posterior,
                posterior_odds=odds,
                genotype_pairs=call.top_class.genotype_pairs,
                assigned_reads=call.assigned_reads,
                unresolved_fraction=call.unresolved_fraction,
                execution_backend=backend_used,
            )
        )
    output.sort(key=lambda item: tuple(item.parameters[name] for name in sorted(item.parameters)))
    return output


def stability_summary(
    scenarios: Sequence[StabilityScenario],
    *,
    base_config: GenotypeConfig,
    grid: StabilityGrid,
) -> dict[str, Any]:
    if not scenarios:
        raise StabilityError("cannot summarize an empty stability analysis")
    class_counts = Counter(item.top_class_id for item in scenarios)
    call_counts = Counter(item.call_state for item in scenarios)
    modal_class, modal_count = sorted(
        class_counts.items(), key=lambda item: (-item[1], item[0])
    )[0]
    posteriors = sorted(item.top_posterior for item in scenarios)
    canonical = {
        "base_config": asdict(base_config),
        "axes": {key: list(value) for key, value in sorted(grid.axes.items())},
        "max_scenarios": grid.max_scenarios,
    }
    grid_id = "NGRID_" + sha256(
        json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:16]
    return {
        "schema_version": 1,
        "grid_id": grid_id,
        "scenario_count": len(scenarios),
        "axes": canonical["axes"],
        "modal_top_class_id": modal_class,
        "modal_top_class_fraction": modal_count / len(scenarios),
        "top_class_counts": dict(sorted(class_counts.items())),
        "call_state_counts": dict(sorted(call_counts.items())),
        "resolved_fraction": (
            call_counts.get("RESOLVED_RESEARCH_CALL", 0) / len(scenarios)
        ),
        "top_posterior_min": posteriors[0],
        "top_posterior_median": median(posteriors),
        "top_posterior_max": posteriors[-1],
        "base_config": asdict(base_config),
    }


def write_stability_tsv(
    path: str | Path, scenarios: Sequence[StabilityScenario]
) -> None:
    if not scenarios:
        raise StabilityError("cannot write an empty stability analysis")
    rows = [item.to_dict() for item in scenarios]
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


__all__ = [
    "StabilityError",
    "StabilityGrid",
    "StabilityScenario",
    "load_stability_grid",
    "run_stability_grid",
    "stability_summary",
    "write_stability_tsv",
]
