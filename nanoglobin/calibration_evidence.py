"""Truth-scoped control evidence for NanoGlobin model calibration."""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Sequence

from .calibration_types import (
    _ELIGIBLE_ASSAY_SCOPES,
    CalibrationError,
    CalibrationSettings,
    ControlEvidence,
    ControlManifestRow,
)
from .genotype import (
    CompiledGenotypeSpace,
    GenotypeConfig,
    MoleculeEvidence,
    read_molecule_evidence,
)


def _weighted_evidence(
    evidence: Sequence[MoleculeEvidence], config: GenotypeConfig
) -> tuple[list[tuple[MoleculeEvidence, float]], float]:
    weighted: list[tuple[MoleculeEvidence, float]] = []
    total = 0.0
    for item in evidence:
        if item.status != "complete" or not item.assigned_product_id:
            continue
        weight = float(config.assignment_weights.get(item.assignment_state, 0.0))
        if weight <= 0:
            continue
        weighted.append((item, weight))
        total += weight
    if total <= 0:
        return [], 0.0
    scale = min(1.0, float(config.effective_read_cap) / total)
    return [(item, weight * scale) for item, weight in weighted], total * scale


def _expected_for_pair(
    space: CompiledGenotypeSpace, pair: tuple[str, str]
) -> tuple[Counter[str], set[tuple[str, str]], set[str]]:
    unknown = sorted(set(pair) - set(space.haplotype_ids))
    if unknown:
        raise CalibrationError(f"unknown calibration haplotype(s): {', '.join(unknown)}")
    multiplicity: Counter[str] = Counter()
    exact_keys: set[tuple[str, str]] = set()
    required: set[str] = set()
    for haplotype_id in pair:
        for record in space.products_by_haplotype[haplotype_id]:
            multiplicity[record.key.product_id] += 1
            exact_keys.add((record.key.product_id, record.key.sequence_sha256))
            if record.required:
                required.add(record.key.product_id)
    if not multiplicity:
        raise CalibrationError(
            f"calibration genotype {pair[0]}/{pair[1]} produces no compiled products"
        )
    return multiplicity, exact_keys, required


def build_control_evidence(
    space: CompiledGenotypeSpace,
    rows: Sequence[ControlManifestRow],
    base_config: GenotypeConfig,
    settings: CalibrationSettings,
) -> list[ControlEvidence]:
    settings.validate()
    output: list[ControlEvidence] = []
    for row in rows:
        pair = row.genotype_pair
        multiplicity, exact_keys, required = _expected_for_pair(space, pair)
        evidence = read_molecule_evidence(row.molecules_tsv)
        weighted, effective_total = _weighted_evidence(evidence, base_config)
        observed: Counter[str] = Counter()
        unsupported = 0.0
        pair_set = set(pair)
        for item, weight in weighted:
            key = (item.assigned_product_id, item.best_sequence_sha256)
            if item.best_sequence_sha256:
                supported = key in exact_keys
            else:
                supported = (
                    item.assigned_product_id in multiplicity
                    and bool(pair_set.intersection(item.candidate_haplotypes))
                )
            if supported:
                observed[item.assigned_product_id] += weight
            else:
                unsupported += weight
        supported_total = sum(observed.values())
        exclusion_reason = ""
        comparator_enabled = (
            settings.allow_comparator_controls
            and row.truth_status == "comparator_only"
            and row.assay_scope == "caller_comparator"
        )
        if (
            row.truth_status == "comparator_only"
            and not settings.allow_comparator_controls
        ):
            exclusion_reason = "comparator_only control is disabled by default"
        elif (
            row.assay_scope not in _ELIGIBLE_ASSAY_SCOPES
            and not comparator_enabled
        ):
            exclusion_reason = (
                f"assay_scope={row.assay_scope!r} is not assay calibration evidence"
            )
        elif supported_total < settings.min_effective_reads:
            exclusion_reason = (
                f"supported effective reads {supported_total:.3f} are below "
                f"{settings.min_effective_reads:.3f}"
            )
        output.append(
            ControlEvidence(
                sample_id=row.sample_id,
                stratum=row.stratum,
                genotype_pair=pair,
                truth_status=row.truth_status,
                truth_source=row.truth_source,
                assay_scope=row.assay_scope,
                expected_multiplicity=dict(multiplicity),
                required_products=frozenset(required),
                observed_counts=dict(observed),
                total_effective_weight=effective_total,
                supported_effective_weight=supported_total,
                unsupported_effective_weight=unsupported,
                eligible=not exclusion_reason,
                exclusion_reason=exclusion_reason,
            )
        )
    return output


def _connected_components(
    controls: Sequence[ControlEvidence],
) -> list[tuple[str, ...]]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for control in controls:
        if not control.eligible:
            continue
        products = sorted(control.expected_multiplicity)
        for product_id in products:
            adjacency[product_id].add(product_id)
        for left in products:
            adjacency[left].update(products)
    remaining = set(adjacency)
    components: list[tuple[str, ...]] = []
    while remaining:
        start = min(remaining)
        stack = [start]
        seen: set[str] = set()
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            stack.extend(sorted(adjacency[node] - seen, reverse=True))
        remaining -= seen
        components.append(tuple(sorted(seen)))
    return sorted(components)
