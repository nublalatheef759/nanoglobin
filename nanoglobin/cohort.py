"""Provenance-preserving normalisation of report-derived cohort tables.

Clinical report exports are selective observations, not complete variant callsets.
This module keeps that distinction explicit while splitting a wide spreadsheet
into test episodes, reported findings, and dated phenotype measurements.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import csv
import re


class CohortContractError(ValueError):
    """Raised when the cohort column contract is invalid or cannot be applied."""


@dataclass(frozen=True)
class LocusColumns:
    locus: str
    result_column: str
    zygosity_column: str | None = None
    details_column: str | None = None


@dataclass(frozen=True)
class PhenotypeColumn:
    measurement: str
    column: str
    unit: str = ""
    date_column: str | None = None


@dataclass(frozen=True)
class CohortContract:
    contract_id: str
    person_id_column: str | None
    episode_id_column: str | None
    shared_columns: Mapping[str, str]
    loci: tuple[LocusColumns, ...]
    phenotypes: tuple[PhenotypeColumn, ...]
    not_tested_values: frozenset[str]
    missing_values: frozenset[str]


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _normalised_token(value: Any) -> str:
    return re.sub(r"[\s_-]+", "", _clean_text(value).casefold())


def cohort_contract_from_mapping(raw: Mapping[str, Any]) -> CohortContract:
    if not isinstance(raw, Mapping):
        raise CohortContractError("cohort contract must be a mapping")
    contract_id = _clean_text(raw.get("contract_id"))
    if not contract_id:
        raise CohortContractError("contract_id is required")

    ids = raw.get("identifiers") or {}
    if not isinstance(ids, Mapping):
        raise CohortContractError("identifiers must be a mapping")
    person_id_column = _clean_text(ids.get("person_id")) or None
    episode_id_column = _clean_text(ids.get("episode_id")) or None

    shared = raw.get("shared_columns") or {}
    if not isinstance(shared, Mapping):
        raise CohortContractError("shared_columns must be a mapping")
    shared_columns = {
        _clean_text(name): _clean_text(column)
        for name, column in shared.items()
        if _clean_text(name) and _clean_text(column)
    }

    loci_raw = raw.get("loci")
    if not isinstance(loci_raw, Mapping) or not loci_raw:
        raise CohortContractError("loci must be a non-empty mapping")
    loci: list[LocusColumns] = []
    for locus, spec_any in loci_raw.items():
        locus_name = _clean_text(locus).upper()
        if not locus_name:
            raise CohortContractError("locus names cannot be empty")
        if not isinstance(spec_any, Mapping):
            raise CohortContractError(f"locus {locus_name} must be a mapping")
        result_column = _clean_text(spec_any.get("result"))
        if not result_column:
            raise CohortContractError(f"locus {locus_name} requires a result column")
        loci.append(
            LocusColumns(
                locus=locus_name,
                result_column=result_column,
                zygosity_column=_clean_text(spec_any.get("zygosity")) or None,
                details_column=_clean_text(spec_any.get("details")) or None,
            )
        )

    phenotypes_raw = raw.get("phenotypes") or {}
    if not isinstance(phenotypes_raw, Mapping):
        raise CohortContractError("phenotypes must be a mapping")
    phenotypes: list[PhenotypeColumn] = []
    for measurement, spec_any in phenotypes_raw.items():
        name = _clean_text(measurement)
        if isinstance(spec_any, str):
            column, unit, date_column = _clean_text(spec_any), "", None
        elif isinstance(spec_any, Mapping):
            column = _clean_text(spec_any.get("column"))
            unit = _clean_text(spec_any.get("unit"))
            date_column = _clean_text(spec_any.get("date_column")) or None
        else:
            raise CohortContractError(f"phenotype {name} must be a string or mapping")
        if not name or not column:
            raise CohortContractError("phenotype names and columns cannot be empty")
        phenotypes.append(PhenotypeColumn(name, column, unit, date_column))

    missing = raw.get("missing_values", ["", "NA", "N/A", "NULL", "NONE"])
    not_tested = raw.get("not_tested_values", ["not tested", "not_tested"])
    if isinstance(missing, str) or not isinstance(missing, Sequence):
        raise CohortContractError("missing_values must be a list")
    if isinstance(not_tested, str) or not isinstance(not_tested, Sequence):
        raise CohortContractError("not_tested_values must be a list")

    return CohortContract(
        contract_id=contract_id,
        person_id_column=person_id_column,
        episode_id_column=episode_id_column,
        shared_columns=shared_columns,
        loci=tuple(loci),
        phenotypes=tuple(phenotypes),
        missing_values=frozenset(_normalised_token(value) for value in missing),
        not_tested_values=frozenset(_normalised_token(value) for value in not_tested),
    )


def load_cohort_contract(path: str | Path) -> CohortContract:
    profile_path = Path(path)
    text = profile_path.read_text(encoding="utf-8")
    if profile_path.suffix.lower() == ".json":
        import json

        raw = json.loads(text)
    else:
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("PyYAML is required to load YAML cohort contracts") from exc
        raw = yaml.safe_load(text)
    return cohort_contract_from_mapping(raw)


def required_source_columns(contract: CohortContract) -> set[str]:
    columns = set(contract.shared_columns.values())
    if contract.person_id_column:
        columns.add(contract.person_id_column)
    if contract.episode_id_column:
        columns.add(contract.episode_id_column)
    for locus in contract.loci:
        columns.add(locus.result_column)
        if locus.zygosity_column:
            columns.add(locus.zygosity_column)
        if locus.details_column:
            columns.add(locus.details_column)
    for phenotype in contract.phenotypes:
        columns.add(phenotype.column)
        if phenotype.date_column:
            columns.add(phenotype.date_column)
    return columns


def _is_missing(value: Any, contract: CohortContract) -> bool:
    return _normalised_token(value) in contract.missing_values


def _is_not_tested(value: Any, contract: CohortContract) -> bool:
    return _normalised_token(value) in contract.not_tested_values


def _stable_row_id(source_label: str, row_number: int, prefix: str) -> str:
    digest = sha256(f"{source_label}:{row_number}".encode("utf-8")).hexdigest()[:16]
    return f"{prefix}_{digest}"


def _normalise_age(raw: str) -> tuple[str, str]:
    value = _clean_text(raw)
    if not value:
        return "", ""
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*([ymd]?)", value.casefold())
    if not match:
        return "", ""
    number = float(match.group(1))
    unit = match.group(2) or "y"
    if unit == "y":
        years = number
    elif unit == "m":
        years = number / 12.0
    else:
        years = number / 365.25
    return f"{years:.6g}", unit


def normalise_rows(
    rows: Iterable[Mapping[str, Any]],
    contract: CohortContract,
    *,
    source_label: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    """Normalise source rows into episodes, findings, and phenotypes.

    Missing locus text is emitted as ``not_reported`` rather than ``reference``.
    This is the critical semantic guard for selectively reported clinical tables.
    """

    episodes: list[dict[str, str]] = []
    findings: list[dict[str, str]] = []
    phenotypes: list[dict[str, str]] = []

    for row_number, row in enumerate(rows, start=2):  # header is line 1
        episode_id = ""
        if contract.episode_id_column:
            episode_id = _clean_text(row.get(contract.episode_id_column))
        if not episode_id:
            episode_id = _stable_row_id(source_label, row_number, "episode")

        person_id = ""
        if contract.person_id_column:
            person_id = _clean_text(row.get(contract.person_id_column))

        episode: dict[str, str] = {
            "episode_id": episode_id,
            "person_id": person_id,
            "source_label": source_label,
            "source_row": str(row_number),
            "contract_id": contract.contract_id,
            "person_linkage_status": "provided" if person_id else "unavailable",
        }
        for canonical, source_column in contract.shared_columns.items():
            raw_value = _clean_text(row.get(source_column))
            episode[canonical] = "" if _is_missing(raw_value, contract) else raw_value
        age_raw = episode.get("age", "")
        age_years, age_source_unit = _normalise_age(age_raw)
        episode["age_years"] = age_years
        episode["age_source_unit"] = age_source_unit
        episodes.append(episode)

        for locus in contract.loci:
            result_raw = _clean_text(row.get(locus.result_column))
            zygosity_raw = (
                _clean_text(row.get(locus.zygosity_column))
                if locus.zygosity_column
                else ""
            )
            details_raw = (
                _clean_text(row.get(locus.details_column))
                if locus.details_column
                else ""
            )
            all_missing = all(
                _is_missing(value, contract)
                for value in (result_raw, zygosity_raw, details_raw)
            )
            if _is_not_tested(result_raw, contract):
                status = "not_tested"
            elif all_missing:
                status = "not_reported"
            else:
                status = "reported"

            findings.append(
                {
                    "finding_id": f"{episode_id}:{locus.locus}",
                    "episode_id": episode_id,
                    "person_id": person_id,
                    "locus": locus.locus,
                    "reporting_status": status,
                    "raw_result": "" if _is_missing(result_raw, contract) else result_raw,
                    "raw_zygosity": ""
                    if _is_missing(zygosity_raw, contract)
                    else zygosity_raw,
                    "raw_details": ""
                    if _is_missing(details_raw, contract)
                    else details_raw,
                    "callset_completeness": "selective_report_only",
                    "reference_genotype_inferred": "false",
                }
            )

        for phenotype in contract.phenotypes:
            raw_value = _clean_text(row.get(phenotype.column))
            if _is_missing(raw_value, contract):
                continue
            measurement_date = (
                _clean_text(row.get(phenotype.date_column))
                if phenotype.date_column
                else ""
            )
            phenotypes.append(
                {
                    "measurement_id": f"{episode_id}:{phenotype.measurement}",
                    "episode_id": episode_id,
                    "person_id": person_id,
                    "measurement": phenotype.measurement,
                    "raw_value": raw_value,
                    "unit": phenotype.unit,
                    "measurement_date": measurement_date,
                    "temporal_context": "unresolved"
                    if not measurement_date
                    else "date_provided_not_adjudicated",
                }
            )

    return episodes, findings, phenotypes


def read_source_csv(
    path: str | Path,
    contract: CohortContract,
) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise CohortContractError(f"{path} has no CSV header")
        missing = sorted(required_source_columns(contract) - set(reader.fieldnames))
        if missing:
            raise CohortContractError(
                f"{path} is missing contract column(s): {', '.join(missing)}"
            )
        return [dict(row) for row in reader]


def write_records(path: str | Path, records: Sequence[Mapping[str, Any]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        output.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    seen: set[str] = set()
    for record in records:
        for field in record:
            if field not in seen:
                fieldnames.append(field)
                seen.add(field)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field, "") for field in fieldnames})
