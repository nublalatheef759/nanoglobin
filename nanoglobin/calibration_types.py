"""Data contracts for NanoGlobin amplicon-model calibration."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from .genotype import GenotypeConfig


class CalibrationError(ValueError):
    """Raised when calibration inputs or fitted parameters are invalid."""


_ELIGIBLE_ASSAY_SCOPES = {"assay_matched", "synthetic_assay"}
_KNOWN_ASSAY_SCOPES = _ELIGIBLE_ASSAY_SCOPES | {
    "wgs",
    "caller_comparator",
    "not_assay_matched",
}
_ALLOWED_TRUTH_STATUS = {
    "certified_reference",
    "orthogonally_resolved",
    "assembly_resolved",
    "synthetic",
    "comparator_only",
}


@dataclass(frozen=True)
class CalibrationSettings:
    """Transparent empirical-Bayes settings for control calibration."""

    min_effective_reads: float = 10.0
    max_iterations: int = 250
    tolerance: float = 1e-9
    efficiency_prior_strength: float = 5.0
    stratum_prior_strength: float = 20.0
    dropout_prior_alpha: float = 1.0
    dropout_prior_beta: float = 19.0
    stratum_dropout_prior_strength: float = 20.0
    artifact_prior_alpha: float = 1.0
    artifact_prior_beta: float = 49.0
    stratum_artifact_prior_strength: float = 50.0
    detection_threshold: float = 0.5
    min_controls_per_stratum: int = 2
    allow_comparator_controls: bool = False

    def validate(self) -> None:
        if self.min_effective_reads <= 0:
            raise CalibrationError("min_effective_reads must be > 0")
        if self.max_iterations < 1:
            raise CalibrationError("max_iterations must be >= 1")
        if self.tolerance <= 0:
            raise CalibrationError("tolerance must be > 0")
        for name in (
            "efficiency_prior_strength",
            "stratum_prior_strength",
            "dropout_prior_alpha",
            "dropout_prior_beta",
            "stratum_dropout_prior_strength",
            "artifact_prior_alpha",
            "artifact_prior_beta",
            "stratum_artifact_prior_strength",
        ):
            if float(getattr(self, name)) <= 0:
                raise CalibrationError(f"{name} must be > 0")
        if self.detection_threshold < 0:
            raise CalibrationError("detection_threshold must be >= 0")
        if self.min_controls_per_stratum < 1:
            raise CalibrationError("min_controls_per_stratum must be >= 1")


@dataclass(frozen=True)
class ControlManifestRow:
    sample_id: str
    molecules_tsv: Path
    haplotype_1: str
    haplotype_2: str
    stratum: str
    truth_status: str
    truth_source: str
    assay_scope: str

    @property
    def genotype_pair(self) -> tuple[str, str]:
        return tuple(sorted((self.haplotype_1, self.haplotype_2)))  # type: ignore[return-value]


@dataclass(frozen=True)
class ControlEvidence:
    sample_id: str
    stratum: str
    genotype_pair: tuple[str, str]
    truth_status: str
    truth_source: str
    assay_scope: str
    expected_multiplicity: Mapping[str, int]
    required_products: frozenset[str]
    observed_counts: Mapping[str, float]
    total_effective_weight: float
    supported_effective_weight: float
    unsupported_effective_weight: float
    eligible: bool
    exclusion_reason: str = ""

    @property
    def detected_products(self) -> frozenset[str]:
        return frozenset(
            product_id
            for product_id, count in self.observed_counts.items()
            if count > 0
        )


@dataclass(frozen=True)
class ComponentFit:
    component_id: str
    products: tuple[str, ...]
    efficiencies: Mapping[str, float]
    observed_weight: Mapping[str, float]
    fitted_exposure: Mapping[str, float]
    controls_expected: Mapping[str, int]
    controls_observed: Mapping[str, int]
    converged: bool
    iterations: int
    identifiability: str


@dataclass(frozen=True)
class CalibrationResult:
    calibration_id: str
    calibration_status: str
    assay_id: str
    assay_version: str
    base_config: GenotypeConfig
    calibrated_config: GenotypeConfig
    settings: CalibrationSettings
    controls: tuple[ControlEvidence, ...]
    components: tuple[ComponentFit, ...]
    global_dropout: Mapping[str, float]
    global_artifact: Mapping[str, float]
    product_dropout: Mapping[str, Mapping[str, float]]
    strata: Mapping[str, Mapping[str, Any]]
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        eligible = [item for item in self.controls if item.eligible]
        excluded = [item for item in self.controls if not item.eligible]
        truth_counts = Counter(item.truth_status for item in eligible)
        scope_counts = Counter(item.assay_scope for item in eligible)
        return {
            "schema_version": 1,
            "calibration_id": self.calibration_id,
            "calibration_status": self.calibration_status,
            "assay_id": self.assay_id,
            "assay_version": self.assay_version,
            "control_summary": {
                "total": len(self.controls),
                "eligible": len(eligible),
                "excluded": len(excluded),
                "truth_status_counts": dict(sorted(truth_counts.items())),
                "assay_scope_counts": dict(sorted(scope_counts.items())),
                "excluded_controls": [
                    {
                        "sample_id": item.sample_id,
                        "reason": item.exclusion_reason,
                        "truth_status": item.truth_status,
                        "assay_scope": item.assay_scope,
                    }
                    for item in excluded
                ],
            },
            "settings": asdict(self.settings),
            "base_model_config": asdict(self.base_config),
            "calibrated_model_config": asdict(self.calibrated_config),
            "parameter_estimates": {
                "dropout": dict(self.global_dropout),
                "unsupported_molecule_mass": dict(self.global_artifact),
            },
            "components": [
                {
                    "component_id": item.component_id,
                    "products": list(item.products),
                    "identifiability": item.identifiability,
                    "converged": item.converged,
                    "iterations": item.iterations,
                    "efficiencies": dict(sorted(item.efficiencies.items())),
                    "observed_weight": dict(sorted(item.observed_weight.items())),
                    "fitted_exposure": dict(sorted(item.fitted_exposure.items())),
                    "controls_expected": dict(sorted(item.controls_expected.items())),
                    "controls_observed": dict(sorted(item.controls_observed.items())),
                }
                for item in self.components
            ],
            "product_dropout": {
                key: dict(value) for key, value in sorted(self.product_dropout.items())
            },
            "strata": {key: dict(value) for key, value in sorted(self.strata.items())},
            "warnings": list(self.warnings),
        }
