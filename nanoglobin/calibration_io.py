"""Calibration manifests, profiles, provenance and tabular outputs."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import csv
import json
import re
from typing import Any, Mapping, Sequence

from .calibration_types import (
    _KNOWN_ASSAY_SCOPES,
    _ALLOWED_TRUTH_STATUS,
    CalibrationError,
    CalibrationResult,
    CalibrationSettings,
    ControlManifestRow,
)
from .genotype import CompiledGenotypeSpace, GenotypeConfig


def load_calibration_settings(path: str | Path | None) -> CalibrationSettings:
    if path is None:
        settings = CalibrationSettings()
        settings.validate()
        return settings
    profile = Path(path)
    text = profile.read_text(encoding="utf-8")
    if profile.suffix.lower() == ".json":
        raw = json.loads(text)
    else:
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("PyYAML is required to load calibration settings") from exc
        raw = yaml.safe_load(text)
    if raw is None:
        raw = {}
    if not isinstance(raw, Mapping):
        raise CalibrationError("calibration settings must be a mapping")
    known = set(CalibrationSettings.__dataclass_fields__)
    unknown = sorted(set(raw) - known)
    if unknown:
        raise CalibrationError(
            f"unknown calibration setting(s): {', '.join(unknown)}"
        )
    settings = CalibrationSettings(**dict(raw))
    settings.validate()
    return settings


def read_control_manifest(path: str | Path) -> list[ControlManifestRow]:
    manifest = Path(path)
    required = {
        "sample_id",
        "molecules_tsv",
        "haplotype_1",
        "haplotype_2",
        "truth_status",
        "truth_source",
        "assay_scope",
    }
    output: list[ControlManifestRow] = []
    seen: set[str] = set()
    with manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise CalibrationError(
                f"{manifest} is missing control column(s): {', '.join(sorted(missing))}"
            )
        for line_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                raise CalibrationError(f"empty sample_id in {manifest}:{line_number}")
            if sample_id in seen:
                raise CalibrationError(f"duplicate sample_id {sample_id!r} in {manifest}")
            seen.add(sample_id)
            truth_status = (row.get("truth_status") or "").strip()
            if truth_status not in _ALLOWED_TRUTH_STATUS:
                raise CalibrationError(
                    f"invalid truth_status {truth_status!r} in {manifest}:{line_number}"
                )
            assay_scope = (row.get("assay_scope") or "").strip()
            if assay_scope not in _KNOWN_ASSAY_SCOPES:
                raise CalibrationError(
                    f"invalid assay_scope {assay_scope!r} in "
                    f"{manifest}:{line_number}"
                )
            molecules_text = (row.get("molecules_tsv") or "").strip()
            if not molecules_text:
                raise CalibrationError(
                    f"empty molecules_tsv in {manifest}:{line_number}"
                )
            molecules = Path(molecules_text)
            if not molecules.is_absolute():
                molecules = manifest.parent / molecules
            haplotype_1 = (row.get("haplotype_1") or "").strip()
            haplotype_2 = (row.get("haplotype_2") or "").strip()
            if not haplotype_1 or not haplotype_2:
                raise CalibrationError(
                    f"empty haplotype in {manifest}:{line_number}"
                )
            truth_source = (row.get("truth_source") or "").strip()
            if not truth_source:
                raise CalibrationError(
                    f"empty truth_source in {manifest}:{line_number}"
                )
            output.append(
                ControlManifestRow(
                    sample_id=sample_id,
                    molecules_tsv=molecules,
                    haplotype_1=haplotype_1,
                    haplotype_2=haplotype_2,
                    stratum=(row.get("stratum") or "global").strip() or "global",
                    truth_status=truth_status,
                    truth_source=truth_source,
                    assay_scope=assay_scope,
                )
            )
    if not output:
        raise CalibrationError(f"control manifest {manifest} contains no controls")
    return output


def product_rows(result: CalibrationResult) -> list[dict[str, Any]]:
    dropout = result.product_dropout
    rows: list[dict[str, Any]] = []
    for component in result.components:
        for product_id in component.products:
            product_dropout = dropout.get(product_id, {})
            rows.append(
                {
                    "component_id": component.component_id,
                    "product_id": product_id,
                    "identifiability": component.identifiability,
                    "relative_efficiency": component.efficiencies[product_id],
                    "observed_weight": component.observed_weight[product_id],
                    "fitted_exposure": component.fitted_exposure[product_id],
                    "controls_expected": component.controls_expected[product_id],
                    "controls_observed": component.controls_observed[product_id],
                    "dropout_events": product_dropout.get("dropout_events", ""),
                    "dropout_trials": product_dropout.get("dropout_trials", ""),
                    "dropout_posterior_mean": product_dropout.get("posterior_mean", ""),
                }
            )
    return rows


def sample_product_rows(result: CalibrationResult) -> list[dict[str, Any]]:
    efficiencies = result.calibrated_config.product_efficiencies
    rows: list[dict[str, Any]] = []
    for control in result.controls:
        if not control.eligible:
            rows.append(
                {
                    "sample_id": control.sample_id,
                    "stratum": control.stratum,
                    "truth_status": control.truth_status,
                    "assay_scope": control.assay_scope,
                    "eligible": False,
                    "exclusion_reason": control.exclusion_reason,
                    "product_id": "",
                    "multiplicity": "",
                    "observed_weight": "",
                    "fitted_weight": "",
                    "residual": "",
                    "required": "",
                    "detected": "",
                }
            )
            continue
        denominator = sum(
            float(multiplicity) * float(efficiencies.get(product_id, 1.0))
            for product_id, multiplicity in control.expected_multiplicity.items()
        )
        for product_id in sorted(control.expected_multiplicity):
            multiplicity = int(control.expected_multiplicity[product_id])
            fitted = (
                control.supported_effective_weight
                * multiplicity
                * float(efficiencies.get(product_id, 1.0))
                / denominator
                if denominator > 0
                else 0.0
            )
            observed = float(control.observed_counts.get(product_id, 0.0))
            rows.append(
                {
                    "sample_id": control.sample_id,
                    "stratum": control.stratum,
                    "truth_status": control.truth_status,
                    "assay_scope": control.assay_scope,
                    "eligible": True,
                    "exclusion_reason": "",
                    "product_id": product_id,
                    "multiplicity": multiplicity,
                    "observed_weight": observed,
                    "fitted_weight": fitted,
                    "residual": observed - fitted,
                    "required": product_id in control.required_products,
                    "detected": observed >= result.settings.detection_threshold,
                }
            )
    return rows


def write_tsv(path: str | Path, rows: Sequence[Mapping[str, Any]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise CalibrationError(f"cannot write empty table to {output}")
    fields = list(rows[0])
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_model_config(path: str | Path, config: GenotypeConfig) -> None:
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("PyYAML is required to write model configuration") from exc
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        yaml.safe_dump(asdict(config), sort_keys=True, default_flow_style=False),
        encoding="utf-8",
    )


def write_stratum_configs(
    directory: str | Path,
    result: CalibrationResult,
) -> list[Path]:
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("PyYAML is required to write stratum configs") from exc
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    base = asdict(result.calibrated_config)
    for stratum, payload in sorted(result.strata.items()):
        if payload.get("status") != "empirical_bayes_stratum_fit":
            continue
        values = dict(base)
        values["product_efficiencies"] = payload["product_efficiencies"]
        values["dropout_probability"] = payload["dropout_probability"]
        values["artifact_probability"] = payload["artifact_probability"]
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", stratum).strip("._") or "stratum"
        path = root / f"{safe_name}.yml"
        path.write_text(
            yaml.safe_dump(values, sort_keys=True, default_flow_style=False),
            encoding="utf-8",
        )
        written.append(path)
    return written


def load_calibration_provenance(
    path: str | Path | None,
    space: CompiledGenotypeSpace,
    config: GenotypeConfig | None = None,
) -> Mapping[str, Any] | None:
    """Load and validate calibration metadata for downstream reporting.

    When ``config`` is supplied, the calibrated configuration embedded in the
    provenance document must exactly match the configuration used for inference.
    This prevents a call from being labelled with calibration produced for a
    different parameter set.
    """

    if path is None:
        return None
    calibration_path = Path(path)
    try:
        payload = json.loads(calibration_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CalibrationError(
            f"cannot read calibration provenance {calibration_path}: {exc}"
        ) from exc
    if not isinstance(payload, Mapping):
        raise CalibrationError("calibration JSON must contain an object")
    if payload.get("assay_id") != space.assay_id:
        raise CalibrationError("calibration assay_id does not match compiled assay")
    if payload.get("assay_version") != space.assay_version:
        raise CalibrationError(
            "calibration assay_version does not match compiled assay"
        )
    calibration_id = payload.get("calibration_id")
    calibration_status = payload.get("calibration_status")
    if not isinstance(calibration_id, str) or not calibration_id:
        raise CalibrationError("calibration JSON lacks calibration_id")
    if not isinstance(calibration_status, str) or not calibration_status:
        raise CalibrationError("calibration JSON lacks calibration_status")
    if config is not None:
        calibrated = payload.get("calibrated_model_config")
        if calibrated != asdict(config):
            raise CalibrationError(
                "calibration provenance does not match the inference model config"
            )
    return {
        "calibration_id": calibration_id,
        "calibration_status": calibration_status,
        "control_summary": payload.get("control_summary"),
        "warnings": payload.get("warnings", []),
    }
