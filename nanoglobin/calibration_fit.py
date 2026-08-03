"""Relative efficiency, dropout and artifact calibration kernels."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict
from hashlib import sha256
import json
from math import exp, log
from typing import Any, Mapping, Sequence

from .calibration_evidence import build_control_evidence, _connected_components
from .calibration_types import (
    CalibrationError,
    CalibrationResult,
    CalibrationSettings,
    ComponentFit,
    ControlEvidence,
    ControlManifestRow,
)
from .genotype import CompiledGenotypeSpace, GenotypeConfig


def _normalize_efficiencies(
    values: Mapping[str, float],
    products: Sequence[str],
    anchor: Mapping[str, float] | None,
) -> dict[str, float]:
    if not products:
        return {}
    if anchor is None:
        logs = [log(max(values[product_id], 1e-12)) for product_id in products]
        shift = exp(sum(logs) / len(logs))
    else:
        logs = [
            log(max(values[product_id], 1e-12) / max(anchor[product_id], 1e-12))
            for product_id in products
        ]
        shift = exp(sum(logs) / len(logs))
    return {
        product_id: max(1e-12, float(values[product_id]) / shift)
        for product_id in products
    }


def _fit_component(
    *,
    component_id: str,
    products: tuple[str, ...],
    controls: Sequence[ControlEvidence],
    prior_efficiencies: Mapping[str, float],
    prior_strength: float,
    max_iterations: int,
    tolerance: float,
    anchor_to_prior: bool,
) -> ComponentFit:
    relevant = [
        control
        for control in controls
        if control.eligible
        and any(control.expected_multiplicity.get(product_id, 0) > 0 for product_id in products)
    ]
    observed_weight = {
        product_id: sum(
            float(control.observed_counts.get(product_id, 0.0)) for control in relevant
        )
        for product_id in products
    }
    controls_expected = {
        product_id: sum(
            control.expected_multiplicity.get(product_id, 0) > 0 for control in relevant
        )
        for product_id in products
    }
    controls_observed = {
        product_id: sum(
            control.observed_counts.get(product_id, 0.0) > 0 for control in relevant
        )
        for product_id in products
    }
    if not relevant:
        return ComponentFit(
            component_id=component_id,
            products=products,
            efficiencies={
                product_id: float(prior_efficiencies.get(product_id, 1.0))
                for product_id in products
            },
            observed_weight=observed_weight,
            fitted_exposure={product_id: 0.0 for product_id in products},
            controls_expected=controls_expected,
            controls_observed=controls_observed,
            converged=True,
            iterations=0,
            identifiability="no_control_overlap_retained_prior",
        )
    if len(products) == 1:
        product_id = products[0]
        value = float(prior_efficiencies.get(product_id, 1.0))
        return ComponentFit(
            component_id=component_id,
            products=products,
            efficiencies={product_id: value},
            observed_weight=observed_weight,
            fitted_exposure={product_id: 0.0},
            controls_expected=controls_expected,
            controls_observed=controls_observed,
            converged=True,
            iterations=0,
            identifiability="single_product_relative_efficiency_unidentifiable",
        )

    prior = {
        product_id: max(1e-12, float(prior_efficiencies.get(product_id, 1.0)))
        for product_id in products
    }
    efficiencies = _normalize_efficiencies(
        prior,
        products,
        prior if anchor_to_prior else None,
    )
    converged = False
    exposures = {product_id: 0.0 for product_id in products}
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        depths: dict[str, float] = {}
        for control in relevant:
            observed_total = sum(
                float(control.observed_counts.get(product_id, 0.0))
                for product_id in products
            )
            denominator = sum(
                float(control.expected_multiplicity.get(product_id, 0))
                * efficiencies[product_id]
                for product_id in products
            )
            if observed_total > 0 and denominator > 0:
                depths[control.sample_id] = observed_total / denominator
        candidate: dict[str, float] = {}
        exposures = {}
        for product_id in products:
            exposure = sum(
                depths.get(control.sample_id, 0.0)
                * float(control.expected_multiplicity.get(product_id, 0))
                for control in relevant
            )
            exposures[product_id] = exposure
            numerator = observed_weight[product_id] + prior_strength * prior[product_id]
            denominator = exposure + prior_strength
            candidate[product_id] = numerator / denominator
        candidate = _normalize_efficiencies(
            candidate,
            products,
            prior if anchor_to_prior else None,
        )
        change = max(
            abs(log(max(candidate[product_id], 1e-12) / efficiencies[product_id]))
            for product_id in products
        )
        efficiencies = candidate
        if change <= tolerance:
            converged = True
            break

    identifiability = "relative_within_connected_component"
    if any(controls_observed[product_id] == 0 for product_id in products):
        identifiability = "weak_relative_fit_unobserved_expected_product"
    return ComponentFit(
        component_id=component_id,
        products=products,
        efficiencies=efficiencies,
        observed_weight=observed_weight,
        fitted_exposure=exposures,
        controls_expected=controls_expected,
        controls_observed=controls_observed,
        converged=converged,
        iterations=iterations,
        identifiability=identifiability,
    )


def _beta_mean(alpha: float, beta: float, events: float, trials: float) -> float:
    if trials < events or events < 0:
        raise CalibrationError("invalid beta-binomial event/trial totals")
    return (alpha + events) / (alpha + beta + trials)


def _dropout_summary(
    controls: Sequence[ControlEvidence],
    *,
    alpha: float,
    beta: float,
    detection_threshold: float,
) -> tuple[
    float | None,
    dict[str, dict[str, float]],
    dict[str, float],
]:
    events = 0.0
    trials = 0.0
    per_product: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0])
    for control in controls:
        if not control.eligible:
            continue
        for product_id in control.required_products:
            if control.expected_multiplicity.get(product_id, 0) <= 0:
                continue
            trials += 1.0
            per_product[product_id][1] += 1.0
            if float(control.observed_counts.get(product_id, 0.0)) < detection_threshold:
                events += 1.0
                per_product[product_id][0] += 1.0
    global_summary = {
        "dropout_events": events,
        "dropout_trials": trials,
        "posterior_alpha": alpha + events,
        "posterior_beta": beta + trials - events,
        "posterior_mean": (
            _beta_mean(alpha, beta, events, trials)
            if trials > 0
            else alpha / (alpha + beta)
        ),
    }
    if trials <= 0:
        return None, {}, global_summary
    output = {
        product_id: {
            "dropout_events": values[0],
            "dropout_trials": values[1],
            "posterior_alpha": alpha + values[0],
            "posterior_beta": beta + values[1] - values[0],
            "posterior_mean": _beta_mean(alpha, beta, values[0], values[1]),
        }
        for product_id, values in sorted(per_product.items())
    }
    return _beta_mean(alpha, beta, events, trials), output, global_summary


def _artifact_summary(
    controls: Sequence[ControlEvidence], *, alpha: float, beta: float
) -> tuple[float | None, dict[str, float]]:
    unsupported = sum(
        control.unsupported_effective_weight for control in controls if control.eligible
    )
    total = sum(
        control.total_effective_weight for control in controls if control.eligible
    )
    summary = {
        "unsupported_effective_weight": unsupported,
        "complete_assigned_effective_weight": total,
        "posterior_alpha": alpha + unsupported,
        "posterior_beta": beta + total - unsupported,
        "posterior_mean": (
            _beta_mean(alpha, beta, unsupported, total)
            if total > 0
            else alpha / (alpha + beta)
        ),
    }
    if total <= 0:
        return None, summary
    return _beta_mean(alpha, beta, unsupported, total), summary


def _calibration_status(controls: Sequence[ControlEvidence]) -> str:
    eligible = [item for item in controls if item.eligible]
    if not eligible:
        return "no_eligible_controls"
    statuses = {item.truth_status for item in eligible}
    scopes = {item.assay_scope for item in eligible}
    if statuses == {"synthetic"} or scopes == {"synthetic_assay"}:
        return "synthetic_fitted_research"
    if "comparator_only" in statuses:
        return "comparator_fitted_research"
    if scopes == {"assay_matched"} and statuses <= {
        "certified_reference",
        "orthogonally_resolved",
        "assembly_resolved",
    }:
        return "assay_matched_control_calibrated_research"
    return "mixed_evidence_calibrated_research"


def calibrate_model(
    space: CompiledGenotypeSpace,
    rows: Sequence[ControlManifestRow],
    *,
    base_config: GenotypeConfig | None = None,
    settings: CalibrationSettings | None = None,
) -> CalibrationResult:
    base_config = base_config or GenotypeConfig()
    base_config.validate()
    settings = settings or CalibrationSettings()
    settings.validate()
    controls = build_control_evidence(space, rows, base_config, settings)
    eligible = [item for item in controls if item.eligible]
    if not eligible:
        reasons = "; ".join(
            f"{item.sample_id}: {item.exclusion_reason}" for item in controls
        )
        raise CalibrationError(f"no eligible calibration controls ({reasons})")

    components_raw = _connected_components(eligible)
    compiled_products = {
        record.key.product_id
        for records in space.products_by_haplotype.values()
        for record in records
    }
    covered_products = {
        product_id
        for component in components_raw
        for product_id in component
    }
    uncovered_products = sorted(compiled_products - covered_products)
    components_raw.extend((product_id,) for product_id in uncovered_products)
    components_raw = sorted(components_raw)

    components: list[ComponentFit] = []
    global_efficiencies = dict(base_config.product_efficiencies)
    warnings: list[str] = []
    if uncovered_products:
        warnings.append(
            "no eligible control genotype covers compiled product(s): "
            + ", ".join(uncovered_products)
        )
    if len(components_raw) > 1:
        warnings.append(
            "relative efficiency scale is not identified between disconnected "
            "product components; each component retains its own declared anchor"
        )
    for index, products in enumerate(components_raw, start=1):
        fit = _fit_component(
            component_id=f"PC{index:03d}",
            products=products,
            controls=eligible,
            prior_efficiencies=base_config.product_efficiencies,
            prior_strength=settings.efficiency_prior_strength,
            max_iterations=settings.max_iterations,
            tolerance=settings.tolerance,
            anchor_to_prior=False,
        )
        components.append(fit)
        global_efficiencies.update(fit.efficiencies)
        if fit.identifiability == "no_control_overlap_retained_prior":
            warnings.append(
                f"{fit.component_id} has no eligible control overlap; retained "
                "declared product-efficiency prior"
            )
        elif len(products) == 1:
            warnings.append(
                f"{fit.component_id} contains only {products[0]!r}; its relative "
                "efficiency is not identifiable from these controls"
            )
        if not fit.converged:
            warnings.append(f"{fit.component_id} did not converge")

    dropout, product_dropout, global_dropout = _dropout_summary(
        eligible,
        alpha=settings.dropout_prior_alpha,
        beta=settings.dropout_prior_beta,
        detection_threshold=settings.detection_threshold,
    )
    if dropout is None:
        dropout = float(base_config.dropout_probability)
        warnings.append(
            "no required-product trials were available; retained base dropout_probability"
        )
    artifact, global_artifact = _artifact_summary(
        eligible,
        alpha=settings.artifact_prior_alpha,
        beta=settings.artifact_prior_beta,
    )
    if artifact is None:
        artifact = float(base_config.artifact_probability)
        warnings.append(
            "no weighted complete molecules were available for artifact calibration"
        )
    artifact = min(0.49, max(1e-6, float(artifact)))
    dropout = min(1.0 - 1e-6, max(1e-6, float(dropout)))

    config_values = asdict(base_config)
    config_values["product_efficiencies"] = dict(sorted(global_efficiencies.items()))
    config_values["dropout_probability"] = dropout
    config_values["artifact_probability"] = artifact
    calibrated_config = GenotypeConfig(**config_values)
    calibrated_config.validate()

    strata: dict[str, dict[str, Any]] = {}
    for stratum in sorted({item.stratum for item in eligible}):
        subset = [item for item in eligible if item.stratum == stratum]
        if len(subset) < settings.min_controls_per_stratum:
            strata[stratum] = {
                "status": "insufficient_controls",
                "control_count": len(subset),
            }
            continue
        stratum_efficiencies = dict(global_efficiencies)
        stratum_components: list[dict[str, Any]] = []
        for fit in components:
            stratum_fit = _fit_component(
                component_id=fit.component_id,
                products=fit.products,
                controls=subset,
                prior_efficiencies=global_efficiencies,
                prior_strength=settings.stratum_prior_strength,
                max_iterations=settings.max_iterations,
                tolerance=settings.tolerance,
                anchor_to_prior=True,
            )
            stratum_efficiencies.update(stratum_fit.efficiencies)
            stratum_components.append(
                {
                    "component_id": stratum_fit.component_id,
                    "efficiencies": dict(sorted(stratum_fit.efficiencies.items())),
                    "identifiability": stratum_fit.identifiability,
                    "converged": stratum_fit.converged,
                    "iterations": stratum_fit.iterations,
                }
            )
        stratum_dropout, _, stratum_dropout_summary = _dropout_summary(
            subset,
            alpha=dropout * settings.stratum_dropout_prior_strength,
            beta=(1.0 - dropout) * settings.stratum_dropout_prior_strength,
            detection_threshold=settings.detection_threshold,
        )
        stratum_artifact, stratum_artifact_summary = _artifact_summary(
            subset,
            alpha=artifact * settings.stratum_artifact_prior_strength,
            beta=(1.0 - artifact) * settings.stratum_artifact_prior_strength,
        )
        strata[stratum] = {
            "status": "empirical_bayes_stratum_fit",
            "control_count": len(subset),
            "product_efficiencies": dict(sorted(stratum_efficiencies.items())),
            "dropout_probability": (
                dropout if stratum_dropout is None else float(stratum_dropout)
            ),
            "artifact_probability": min(
                0.49,
                max(1e-6, artifact if stratum_artifact is None else float(stratum_artifact)),
            ),
            "dropout_diagnostics": stratum_dropout_summary,
            "artifact_diagnostics": stratum_artifact_summary,
            "components": stratum_components,
        }

    canonical = {
        "assay_id": space.assay_id,
        "assay_version": space.assay_version,
        "eligible_controls": [
            {
                "sample_id": item.sample_id,
                "pair": list(item.genotype_pair),
                "truth_status": item.truth_status,
                "truth_source": item.truth_source,
                "assay_scope": item.assay_scope,
                "stratum": item.stratum,
            }
            for item in eligible
        ],
        "config": asdict(calibrated_config),
        "settings": asdict(settings),
    }
    digest = sha256(
        json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:16]
    return CalibrationResult(
        calibration_id=f"NGCAL_{digest}",
        calibration_status=_calibration_status(controls),
        assay_id=space.assay_id,
        assay_version=space.assay_version,
        base_config=base_config,
        calibrated_config=calibrated_config,
        settings=settings,
        controls=tuple(controls),
        components=tuple(components),
        global_dropout=global_dropout,
        global_artifact=global_artifact,
        product_dropout=product_dropout,
        strata=strata,
        warnings=tuple(warnings),
    )
