"""Candidate chromosome-haplotype posterior for declared amplicon assays.

The candidate space is a pair of chromosome haplotypes. Each chromosome
haplotype may itself encode any HBA copy number, copy order, hybrid gene or linked
small variants. Therefore a variable-copy HBA model does not require changing
human ploidy: it requires sequence-resolved chromosome-haplotype candidates and a
posterior over their diploid pairs.

The implementation combines three transparent terms:

* molecule-to-compiled-product compatibility;
* a Dirichlet-multinomial product-composition score with explicit efficiencies;
* required-product/dropout evidence.

Named genotype pairs that compile to the same count-aware product signature are
collapsed before assigning the default prior. Adding synonymous catalogue labels
therefore cannot manufacture posterior mass. The resulting posterior is a
research posterior conditional on the candidate catalogue and declared model; it
is not claimed to be clinically calibrated before truth-matched validation.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from itertools import combinations_with_replacement
from math import exp, lgamma, log
from pathlib import Path
from typing import Any, Mapping, Sequence
import csv
import hashlib
import json


class GenotypeModelError(ValueError):
    """Raised when compiled products, molecule evidence or config are invalid."""


@dataclass(frozen=True, order=True)
class ProductKey:
    product_id: str
    sequence_sha256: str


@dataclass(frozen=True)
class CompiledProductRecord:
    haplotype_id: str
    key: ProductKey
    required: bool


@dataclass(frozen=True)
class CompiledGenotypeSpace:
    assay_id: str
    assay_version: str
    haplotype_ids: tuple[str, ...]
    products_by_haplotype: Mapping[str, tuple[CompiledProductRecord, ...]]


@dataclass(frozen=True)
class MoleculeEvidence:
    read_id: str
    status: str
    assigned_product_id: str
    assignment_state: str
    best_sequence_sha256: str
    candidate_haplotypes: tuple[str, ...]
    edit_rate: float | None


@dataclass(frozen=True)
class GenotypeConfig:
    """Explicit parameters for the bounded candidate posterior."""

    min_assigned_reads: int = 3
    effective_read_cap: float = 200.0
    effective_count_cap: float = 200.0
    artifact_probability: float = 0.02
    product_background_mass: float = 0.01
    count_concentration: float = 20.0
    minimum_dirichlet_alpha: float = 0.05
    dropout_probability: float = 0.05
    resolved_posterior: float = 0.80
    resolved_odds: float = 3.0
    max_unresolved_fraction: float = 0.80
    assignment_weights: Mapping[str, float] = field(
        default_factory=lambda: {
            "unique_sequence": 1.0,
            "equivalent_haplotypes": 1.0,
            "low_margin": 0.35,
            "ambiguous_product_definition": 0.25,
            "poor_sequence_fit": 0.10,
            "no_compiled_sequence": 0.0,
            "": 0.0,
        }
    )
    product_efficiencies: Mapping[str, float] = field(default_factory=dict)
    haplotype_priors: Mapping[str, float] = field(default_factory=dict)

    def validate(self) -> None:
        if self.min_assigned_reads < 1:
            raise GenotypeModelError("min_assigned_reads must be >= 1")
        if self.effective_read_cap <= 0 or self.effective_count_cap <= 0:
            raise GenotypeModelError("effective evidence caps must be > 0")
        if not 0 < self.artifact_probability < 0.5:
            raise GenotypeModelError("artifact_probability must be in (0, 0.5)")
        if self.product_background_mass <= 0:
            raise GenotypeModelError("product_background_mass must be > 0")
        if self.count_concentration <= 0 or self.minimum_dirichlet_alpha <= 0:
            raise GenotypeModelError("Dirichlet parameters must be > 0")
        if not 0 < self.dropout_probability < 1:
            raise GenotypeModelError("dropout_probability must be in (0, 1)")
        if not 0 < self.resolved_posterior <= 1:
            raise GenotypeModelError("resolved_posterior must be in (0, 1]")
        if self.resolved_odds < 1:
            raise GenotypeModelError("resolved_odds must be >= 1")
        if not 0 <= self.max_unresolved_fraction <= 1:
            raise GenotypeModelError("max_unresolved_fraction must be in [0, 1]")
        for state, weight in self.assignment_weights.items():
            if weight < 0:
                raise GenotypeModelError(
                    f"assignment weight for {state!r} must be >= 0"
                )
        for product_id, efficiency in self.product_efficiencies.items():
            if efficiency <= 0:
                raise GenotypeModelError(
                    f"product efficiency for {product_id!r} must be > 0"
                )
        if self.haplotype_priors and any(
            value <= 0 for value in self.haplotype_priors.values()
        ):
            raise GenotypeModelError("haplotype priors must be > 0")


@dataclass(frozen=True)
class ObservableGenotypeClass:
    class_id: str
    genotype_pairs: tuple[tuple[str, str], ...]
    signature: tuple[tuple[str, str, int], ...]


@dataclass(frozen=True)
class GenotypeClassScore:
    class_id: str
    genotype_pairs: tuple[tuple[str, str], ...]
    signature: tuple[tuple[str, str, int], ...]
    log_read_likelihood: float
    log_count_likelihood: float
    log_dropout_likelihood: float
    log_prior: float
    log_score: float
    posterior: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "class_id": self.class_id,
            "genotype_pairs": [list(pair) for pair in self.genotype_pairs],
            "signature": [list(item) for item in self.signature],
            "log_read_likelihood": self.log_read_likelihood,
            "log_count_likelihood": self.log_count_likelihood,
            "log_dropout_likelihood": self.log_dropout_likelihood,
            "log_prior": self.log_prior,
            "log_score": self.log_score,
            "posterior": self.posterior,
        }


@dataclass(frozen=True)
class GenotypeCall:
    call_state: str
    reason: str
    top_class: GenotypeClassScore | None
    second_class: GenotypeClassScore | None
    assigned_reads: int
    useful_effective_reads: float
    reads_total: int
    unresolved_fraction: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "call_state": self.call_state,
            "reason": self.reason,
            "top_class": self.top_class.to_dict() if self.top_class else None,
            "second_class": self.second_class.to_dict() if self.second_class else None,
            "assigned_reads": self.assigned_reads,
            "useful_effective_reads": self.useful_effective_reads,
            "reads_total": self.reads_total,
            "unresolved_fraction": self.unresolved_fraction,
        }


def load_compiled_space(path: str | Path) -> CompiledGenotypeSpace:
    """Load the JSON emitted by ``compile_assay.py``."""

    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GenotypeModelError(f"cannot read compiled assay {path}: {exc}") from exc
    assay = payload.get("assay")
    haplotypes = payload.get("haplotypes")
    products = payload.get("compiled_products")
    if not isinstance(assay, Mapping):
        raise GenotypeModelError("compiled assay lacks assay metadata")
    if not isinstance(haplotypes, Sequence) or isinstance(haplotypes, (str, bytes)):
        raise GenotypeModelError("compiled assay lacks haplotypes")
    if not isinstance(products, Sequence) or isinstance(products, (str, bytes)):
        raise GenotypeModelError("compiled assay lacks compiled_products")
    assay_id = str(assay.get("assay_id", "")).strip()
    assay_version = str(assay.get("version", "")).strip()
    if not assay_id or not assay_version:
        raise GenotypeModelError("compiled assay lacks assay id/version")

    haplotype_ids: list[str] = []
    for index, raw in enumerate(haplotypes):
        if not isinstance(raw, Mapping):
            raise GenotypeModelError(f"haplotypes[{index}] is not an object")
        haplotype_id = str(raw.get("haplotype_id", "")).strip()
        if not haplotype_id:
            raise GenotypeModelError(f"haplotypes[{index}] lacks haplotype_id")
        haplotype_ids.append(haplotype_id)
    if len(set(haplotype_ids)) != len(haplotype_ids):
        raise GenotypeModelError("compiled assay contains duplicate haplotype IDs")

    products_by_haplotype: dict[str, list[CompiledProductRecord]] = {
        haplotype_id: [] for haplotype_id in haplotype_ids
    }
    for index, raw in enumerate(products):
        if not isinstance(raw, Mapping):
            raise GenotypeModelError(f"compiled_products[{index}] is not an object")
        try:
            haplotype_id = str(raw["haplotype_id"])
            product_id = str(raw["product_id"])
            sequence_hash = str(raw["sequence_sha256"])
        except KeyError as exc:
            raise GenotypeModelError(
                f"compiled_products[{index}] lacks {exc.args[0]}"
            ) from exc
        if haplotype_id not in products_by_haplotype:
            raise GenotypeModelError(
                f"compiled product references unknown haplotype {haplotype_id!r}"
            )
        if not product_id or not sequence_hash:
            raise GenotypeModelError(f"compiled_products[{index}] has an empty key")
        products_by_haplotype[haplotype_id].append(
            CompiledProductRecord(
                haplotype_id=haplotype_id,
                key=ProductKey(product_id, sequence_hash),
                required=bool(raw.get("required", False)),
            )
        )

    return CompiledGenotypeSpace(
        assay_id=assay_id,
        assay_version=assay_version,
        haplotype_ids=tuple(sorted(haplotype_ids)),
        products_by_haplotype={
            haplotype_id: tuple(
                sorted(
                    records,
                    key=lambda item: (item.key.product_id, item.key.sequence_sha256),
                )
            )
            for haplotype_id, records in products_by_haplotype.items()
        },
    )


def read_molecule_evidence(path: str | Path) -> list[MoleculeEvidence]:
    """Read molecule evidence emitted by ``admit_amplicon_reads.py``."""

    required = {
        "read_id",
        "status",
        "assigned_product_id",
        "assignment_state",
        "best_sequence_sha256",
        "candidate_haplotypes",
        "edit_rate",
    }
    output: list[MoleculeEvidence] = []
    seen_reads: set[str] = set()
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise GenotypeModelError(
                f"{path} is missing molecule column(s): {', '.join(sorted(missing))}"
            )
        for line_number, row in enumerate(reader, start=2):
            read_id = (row.get("read_id") or "").strip()
            if not read_id:
                raise GenotypeModelError(f"empty read_id in {path}:{line_number}")
            if read_id in seen_reads:
                raise GenotypeModelError(f"duplicate read_id {read_id!r} in {path}")
            seen_reads.add(read_id)
            edit_text = (row.get("edit_rate") or "").strip()
            try:
                edit_rate = float(edit_text) if edit_text else None
            except ValueError as exc:
                raise GenotypeModelError(
                    f"invalid edit_rate in {path}:{line_number}: {edit_text!r}"
                ) from exc
            output.append(
                MoleculeEvidence(
                    read_id=read_id,
                    status=(row.get("status") or "").strip(),
                    assigned_product_id=(row.get("assigned_product_id") or "").strip(),
                    assignment_state=(row.get("assignment_state") or "").strip(),
                    best_sequence_sha256=(row.get("best_sequence_sha256") or "").strip(),
                    candidate_haplotypes=tuple(
                        sorted(
                            value
                            for value in (row.get("candidate_haplotypes") or "").split(";")
                            if value
                        )
                    ),
                    edit_rate=edit_rate,
                )
            )
    return output


def genotype_signature(
    space: CompiledGenotypeSpace, pair: tuple[str, str]
) -> tuple[tuple[str, str, int], ...]:
    counts: Counter[ProductKey] = Counter()
    for haplotype_id in pair:
        counts.update(record.key for record in space.products_by_haplotype[haplotype_id])
    return tuple(
        sorted(
            (key.product_id, key.sequence_sha256, count)
            for key, count in counts.items()
        )
    )


def observable_classes(space: CompiledGenotypeSpace) -> list[ObservableGenotypeClass]:
    """Collapse named genotype pairs by count-aware product signature."""

    grouped: dict[
        tuple[tuple[str, str, int], ...], list[tuple[str, str]]
    ] = defaultdict(list)
    for pair in combinations_with_replacement(space.haplotype_ids, 2):
        grouped[genotype_signature(space, pair)].append(pair)
    output: list[ObservableGenotypeClass] = []
    for signature, pairs in sorted(grouped.items(), key=lambda item: (item[0], item[1])):
        digest = json.dumps(signature, separators=(",", ":")).encode("utf-8")
        class_id = "OGC_" + hashlib.sha256(digest).hexdigest()[:12]
        output.append(
            ObservableGenotypeClass(
                class_id=class_id,
                genotype_pairs=tuple(sorted(pairs)),
                signature=signature,
            )
        )
    return output


def _logsumexp(values: Sequence[float]) -> float:
    if not values:
        raise GenotypeModelError("cannot normalize an empty score set")
    maximum = max(values)
    return maximum + log(sum(exp(value - maximum) for value in values))


def _state_weight(evidence: MoleculeEvidence, config: GenotypeConfig) -> float:
    if evidence.status != "complete" or not evidence.assigned_product_id:
        return 0.0
    return float(config.assignment_weights.get(evidence.assignment_state, 0.0))


def _effective_weights(
    evidence: Sequence[MoleculeEvidence], config: GenotypeConfig
) -> tuple[list[tuple[MoleculeEvidence, float]], float]:
    weighted = [(item, _state_weight(item, config)) for item in evidence]
    total = sum(weight for _, weight in weighted)
    if total <= 0:
        return [], 0.0
    scale = min(1.0, config.effective_read_cap / total)
    return [(item, weight * scale) for item, weight in weighted if weight > 0], total * scale


def _predicted_key_set(
    signature: tuple[tuple[str, str, int], ...]
) -> set[ProductKey]:
    return {
        ProductKey(product_id, sequence_hash)
        for product_id, sequence_hash, _ in signature
    }


def _predicted_product_ids(
    signature: tuple[tuple[str, str, int], ...]
) -> set[str]:
    return {product_id for product_id, _, _ in signature}


def _read_log_likelihood(
    genotype_class: ObservableGenotypeClass,
    weighted_evidence: Sequence[tuple[MoleculeEvidence, float]],
    config: GenotypeConfig,
) -> float:
    predicted_keys = _predicted_key_set(genotype_class.signature)
    predicted_products = _predicted_product_ids(genotype_class.signature)
    epsilon = config.artifact_probability
    value = 0.0
    for evidence, weight in weighted_evidence:
        key = ProductKey(evidence.assigned_product_id, evidence.best_sequence_sha256)
        pair_haplotypes = {
            haplotype_id
            for pair in genotype_class.genotype_pairs
            for haplotype_id in pair
        }
        candidate_overlap = bool(
            pair_haplotypes.intersection(evidence.candidate_haplotypes)
        )
        if evidence.best_sequence_sha256 and key in predicted_keys:
            probability = 1.0 - epsilon
        elif (
            not evidence.best_sequence_sha256
            and evidence.assigned_product_id in predicted_products
            and candidate_overlap
        ):
            probability = max(epsilon, 0.5 * (1.0 - epsilon))
        elif evidence.assigned_product_id in predicted_products:
            probability = max(epsilon, 0.25 * (1.0 - epsilon))
        else:
            probability = epsilon
        edit_penalty = 0.0
        if evidence.edit_rate is not None:
            edit_penalty = min(1.0, max(0.0, evidence.edit_rate))
        value += weight * (log(max(1e-12, probability)) - edit_penalty)
    return value


def _observed_product_counts(
    weighted_evidence: Sequence[tuple[MoleculeEvidence, float]],
    config: GenotypeConfig,
) -> dict[str, float]:
    counts: Counter[str] = Counter()
    for evidence, weight in weighted_evidence:
        counts[evidence.assigned_product_id] += weight
    total = sum(counts.values())
    if total <= 0:
        return {}
    scale = min(1.0, config.effective_count_cap / total)
    return {product_id: count * scale for product_id, count in counts.items()}


def _predicted_product_multiplicity(
    signature: tuple[tuple[str, str, int], ...]
) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for product_id, _sequence_hash, multiplicity in signature:
        counts[product_id] += multiplicity
    return dict(counts)


def _dirichlet_multinomial_score(
    observed: Mapping[str, float],
    predicted: Mapping[str, int],
    config: GenotypeConfig,
) -> float:
    if not observed:
        return 0.0
    categories = sorted(set(observed) | set(predicted))
    masses: dict[str, float] = {}
    for product_id in categories:
        efficiency = float(config.product_efficiencies.get(product_id, 1.0))
        masses[product_id] = (
            float(predicted.get(product_id, 0)) * efficiency
            + config.product_background_mass
        )
    total_mass = sum(masses.values())
    alpha = {
        product_id: max(
            config.minimum_dirichlet_alpha,
            config.count_concentration * masses[product_id] / total_mass,
        )
        for product_id in categories
    }
    alpha0 = sum(alpha.values())
    n_total = sum(float(observed.get(product_id, 0.0)) for product_id in categories)
    value = lgamma(alpha0) - lgamma(alpha0 + n_total)
    for product_id in categories:
        count = float(observed.get(product_id, 0.0))
        value += lgamma(alpha[product_id] + count) - lgamma(alpha[product_id])
    return value


def _required_product_ids_for_class(
    space: CompiledGenotypeSpace,
    genotype_class: ObservableGenotypeClass,
) -> set[str]:
    required: set[str] = set()
    haplotypes = {
        haplotype_id
        for pair in genotype_class.genotype_pairs
        for haplotype_id in pair
    }
    for haplotype_id in haplotypes:
        for record in space.products_by_haplotype[haplotype_id]:
            if record.required:
                required.add(record.key.product_id)
    return required


def _dropout_log_likelihood(
    space: CompiledGenotypeSpace,
    genotype_class: ObservableGenotypeClass,
    observed: Mapping[str, float],
    config: GenotypeConfig,
) -> float:
    value = 0.0
    for product_id in _required_product_ids_for_class(space, genotype_class):
        value += log(
            1.0 - config.dropout_probability
            if observed.get(product_id, 0.0) > 0
            else config.dropout_probability
        )
    return value


def _class_priors(
    classes: Sequence[ObservableGenotypeClass],
    config: GenotypeConfig,
) -> dict[str, float]:
    if not config.haplotype_priors:
        uniform = 1.0 / len(classes)
        return {item.class_id: uniform for item in classes}

    raw: dict[str, float] = {}
    for item in classes:
        class_mass = 0.0
        for left, right in item.genotype_pairs:
            try:
                left_prior = float(config.haplotype_priors[left])
                right_prior = float(config.haplotype_priors[right])
            except KeyError as exc:
                raise GenotypeModelError(
                    f"missing haplotype prior for {exc.args[0]!r}"
                ) from exc
            class_mass += (
                left_prior * right_prior
                if left == right
                else 2.0 * left_prior * right_prior
            )
        raw[item.class_id] = class_mass
    total = sum(raw.values())
    if total <= 0:
        raise GenotypeModelError("haplotype priors yield zero total mass")
    return {class_id: value / total for class_id, value in raw.items()}


def score_genotypes(
    space: CompiledGenotypeSpace,
    evidence: Sequence[MoleculeEvidence],
    config: GenotypeConfig | None = None,
) -> tuple[list[GenotypeClassScore], float]:
    config = config or GenotypeConfig()
    config.validate()
    classes = observable_classes(space)
    weighted_evidence, useful_effective_reads = _effective_weights(evidence, config)
    observed_counts = _observed_product_counts(weighted_evidence, config)
    priors = _class_priors(classes, config)

    raw: list[tuple[ObservableGenotypeClass, float, float, float, float, float]] = []
    for item in classes:
        read_ll = _read_log_likelihood(item, weighted_evidence, config)
        count_ll = _dirichlet_multinomial_score(
            observed_counts,
            _predicted_product_multiplicity(item.signature),
            config,
        )
        dropout_ll = _dropout_log_likelihood(space, item, observed_counts, config)
        log_prior = log(priors[item.class_id])
        log_score = read_ll + count_ll + dropout_ll + log_prior
        raw.append((item, read_ll, count_ll, dropout_ll, log_prior, log_score))

    normalizer = _logsumexp([item[-1] for item in raw])
    output = [
        GenotypeClassScore(
            class_id=item.class_id,
            genotype_pairs=item.genotype_pairs,
            signature=item.signature,
            log_read_likelihood=read_ll,
            log_count_likelihood=count_ll,
            log_dropout_likelihood=dropout_ll,
            log_prior=log_prior,
            log_score=log_score,
            posterior=exp(log_score - normalizer),
        )
        for item, read_ll, count_ll, dropout_ll, log_prior, log_score in raw
    ]
    output.sort(key=lambda item: (-item.posterior, item.class_id))
    return output, useful_effective_reads


def make_call(
    scores: Sequence[GenotypeClassScore],
    evidence: Sequence[MoleculeEvidence],
    useful_effective_reads: float,
    config: GenotypeConfig | None = None,
) -> GenotypeCall:
    config = config or GenotypeConfig()
    config.validate()
    reads_total = len(evidence)
    assigned_reads = sum(
        item.status == "complete" and bool(item.assigned_product_id)
        for item in evidence
    )
    unresolved_reads = sum(
        item.status != "complete"
        or item.assignment_state
        in {
            "poor_sequence_fit",
            "no_compiled_sequence",
            "ambiguous_product_definition",
        }
        for item in evidence
    )
    unresolved_fraction = unresolved_reads / reads_total if reads_total else 1.0

    top = scores[0] if scores else None
    second = scores[1] if len(scores) > 1 else None
    if assigned_reads < config.min_assigned_reads or useful_effective_reads <= 0:
        return GenotypeCall(
            call_state="NO_CALL_INSUFFICIENT_EVIDENCE",
            reason=(
                f"only {assigned_reads} assigned complete reads; "
                f"minimum is {config.min_assigned_reads}"
            ),
            top_class=top,
            second_class=second,
            assigned_reads=assigned_reads,
            useful_effective_reads=useful_effective_reads,
            reads_total=reads_total,
            unresolved_fraction=unresolved_fraction,
        )
    if unresolved_fraction > config.max_unresolved_fraction:
        return GenotypeCall(
            call_state="NO_CALL_INSUFFICIENT_EVIDENCE",
            reason=(
                f"unresolved molecule fraction {unresolved_fraction:.3f} exceeds "
                f"{config.max_unresolved_fraction:.3f}"
            ),
            top_class=top,
            second_class=second,
            assigned_reads=assigned_reads,
            useful_effective_reads=useful_effective_reads,
            reads_total=reads_total,
            unresolved_fraction=unresolved_fraction,
        )
    if top is None:
        raise GenotypeModelError("no genotype classes were scored")
    odds = top.posterior / max(second.posterior, 1e-300) if second else float("inf")
    if top.posterior < config.resolved_posterior or odds < config.resolved_odds:
        return GenotypeCall(
            call_state="AMBIGUOUS_POSTERIOR",
            reason=(
                f"top posterior={top.posterior:.4f}, odds to second={odds:.3f}; "
                f"required posterior={config.resolved_posterior:.4f} and "
                f"odds={config.resolved_odds:.3f}"
            ),
            top_class=top,
            second_class=second,
            assigned_reads=assigned_reads,
            useful_effective_reads=useful_effective_reads,
            reads_total=reads_total,
            unresolved_fraction=unresolved_fraction,
        )
    if len(top.genotype_pairs) > 1:
        return GenotypeCall(
            call_state="AMBIGUOUS_ASSAY_EQUIVALENT",
            reason=(
                f"{len(top.genotype_pairs)} named genotype pairs share the top "
                "count-aware assay signature"
            ),
            top_class=top,
            second_class=second,
            assigned_reads=assigned_reads,
            useful_effective_reads=useful_effective_reads,
            reads_total=reads_total,
            unresolved_fraction=unresolved_fraction,
        )
    return GenotypeCall(
        call_state="RESOLVED_RESEARCH_CALL",
        reason="one observable genotype class and one named pair exceed thresholds",
        top_class=top,
        second_class=second,
        assigned_reads=assigned_reads,
        useful_effective_reads=useful_effective_reads,
        reads_total=reads_total,
        unresolved_fraction=unresolved_fraction,
    )


def load_config(path: str | Path | None) -> GenotypeConfig:
    if path is None:
        config = GenotypeConfig()
        config.validate()
        return config
    profile_path = Path(path)
    text = profile_path.read_text(encoding="utf-8")
    if profile_path.suffix.lower() == ".json":
        raw = json.loads(text)
    else:
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("PyYAML is required to load genotype config") from exc
        raw = yaml.safe_load(text)
    if raw is None:
        raw = {}
    if not isinstance(raw, Mapping):
        raise GenotypeModelError("genotype config must be a mapping")
    known = {field.name for field in GenotypeConfig.__dataclass_fields__.values()}
    unknown = sorted(set(raw) - known)
    if unknown:
        raise GenotypeModelError(
            f"unknown genotype config field(s): {', '.join(unknown)}"
        )
    config = GenotypeConfig(**dict(raw))
    config.validate()
    return config


def write_scores_tsv(path: str | Path, scores: Sequence[GenotypeClassScore]) -> None:
    fields = [
        "rank",
        "class_id",
        "posterior",
        "genotype_pair_count",
        "genotype_pairs",
        "product_signature",
        "log_read_likelihood",
        "log_count_likelihood",
        "log_dropout_likelihood",
        "log_prior",
        "log_score",
    ]
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rank, item in enumerate(scores, start=1):
            writer.writerow(
                {
                    "rank": rank,
                    "class_id": item.class_id,
                    "posterior": f"{item.posterior:.12g}",
                    "genotype_pair_count": len(item.genotype_pairs),
                    "genotype_pairs": ";".join(
                        f"{left}/{right}" for left, right in item.genotype_pairs
                    ),
                    "product_signature": ";".join(
                        f"{product_id}:{sequence_hash}:{count}"
                        for product_id, sequence_hash, count in item.signature
                    ),
                    "log_read_likelihood": f"{item.log_read_likelihood:.12g}",
                    "log_count_likelihood": f"{item.log_count_likelihood:.12g}",
                    "log_dropout_likelihood": f"{item.log_dropout_likelihood:.12g}",
                    "log_prior": f"{item.log_prior:.12g}",
                    "log_score": f"{item.log_score:.12g}",
                }
            )


def analysis_payload(
    *,
    space: CompiledGenotypeSpace,
    config: GenotypeConfig,
    scores: Sequence[GenotypeClassScore],
    call: GenotypeCall,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "assay_id": space.assay_id,
        "assay_version": space.assay_version,
        "candidate_haplotype_count": len(space.haplotype_ids),
        "observable_genotype_class_count": len(scores),
        "model": {
            "name": "nanoglobin_candidate_product_posterior",
            "calibration_status": "research_uncalibrated",
            "prior_unit": (
                "haplotype-frequency"
                if config.haplotype_priors
                else "observable-signature"
            ),
            "config": asdict(config),
        },
        "call": call.to_dict(),
        "top_classes": [item.to_dict() for item in scores[:10]],
    }


__all__ = [
    "CompiledGenotypeSpace",
    "GenotypeCall",
    "GenotypeClassScore",
    "GenotypeConfig",
    "GenotypeModelError",
    "MoleculeEvidence",
    "ObservableGenotypeClass",
    "analysis_payload",
    "genotype_signature",
    "load_compiled_space",
    "load_config",
    "make_call",
    "observable_classes",
    "read_molecule_evidence",
    "score_genotypes",
    "write_scores_tsv",
]
