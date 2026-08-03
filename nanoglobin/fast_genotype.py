"""Automatic Python/Cython backends for chromosome-haplotype scoring.

RapidFuzz already executes edit-distance work in native code.  This module targets
the remaining potentially quadratic stage: scoring every assay-observable genotype
class against molecule evidence and product-composition terms.

The Cython path first groups equivalent molecule observations.  A 20,000-read PCR
run can therefore collapse to a small number of product/sequence/candidate states
before the class-by-evidence kernel runs.  The reference Python implementation is
retained as a semantic oracle and automatic fallback.
"""

from __future__ import annotations

from collections import defaultdict
from math import exp, log
from typing import Literal, Sequence

import numpy as np

from .genotype import (
    CompiledGenotypeSpace,
    GenotypeClassScore,
    GenotypeConfig,
    GenotypeModelError,
    MoleculeEvidence,
    ObservableGenotypeClass,
    ProductKey,
    _class_priors,
    _effective_weights,
    _logsumexp,
    _observed_product_counts,
    _required_product_ids_for_class,
    observable_classes,
    score_genotypes,
)

BackendName = Literal["auto", "python", "cython"]


def cython_backend_available() -> bool:
    try:
        from . import _genotype_fast  # noqa: F401
    except ImportError:
        return False
    return True


def _fallback_or_raise(
    *,
    backend: BackendName,
    reason: str,
    space: CompiledGenotypeSpace,
    evidence: Sequence[MoleculeEvidence],
    config: GenotypeConfig,
) -> tuple[list[GenotypeClassScore], float, str]:
    if backend == "cython":
        raise GenotypeModelError(f"Cython genotype backend unavailable: {reason}")
    scores, effective = score_genotypes(space, evidence, config)
    return scores, effective, f"python_reference ({reason})"


def _haplotype_mask(
    haplotypes: Sequence[str], index: dict[str, int]
) -> np.uint64:
    value = 0
    for haplotype_id in haplotypes:
        position = index.get(haplotype_id)
        if position is not None:
            value |= 1 << position
    return np.uint64(value)


def _prepare_dense_problem(
    *,
    space: CompiledGenotypeSpace,
    classes: Sequence[ObservableGenotypeClass],
    evidence: Sequence[MoleculeEvidence],
    config: GenotypeConfig,
):
    weighted_evidence, useful_effective_reads = _effective_weights(evidence, config)
    observed = _observed_product_counts(weighted_evidence, config)

    product_ids = sorted(
        set(observed)
        | {
            product_id
            for item in classes
            for product_id, _sequence_hash, _multiplicity in item.signature
        }
    )
    product_index = {product_id: index for index, product_id in enumerate(product_ids)}

    exact_keys = sorted(
        {
            ProductKey(product_id, sequence_hash)
            for item in classes
            for product_id, sequence_hash, _multiplicity in item.signature
        }
    )
    exact_key_index = {key: index for index, key in enumerate(exact_keys)}
    haplotype_ids = tuple(sorted(space.haplotype_ids))
    haplotype_index = {
        haplotype_id: index for index, haplotype_id in enumerate(haplotype_ids)
    }

    class_count = len(classes)
    product_count = len(product_ids)
    key_count = len(exact_keys)
    class_exact = np.zeros((class_count, key_count), dtype=np.uint8)
    class_product = np.zeros((class_count, product_count), dtype=np.uint8)
    class_required = np.zeros((class_count, product_count), dtype=np.uint8)
    class_multiplicity = np.zeros((class_count, product_count), dtype=np.int64)
    class_haplotype_masks = np.zeros(class_count, dtype=np.uint64)

    for class_index, item in enumerate(classes):
        class_haplotypes = sorted(
            {
                haplotype_id
                for pair in item.genotype_pairs
                for haplotype_id in pair
            }
        )
        class_haplotype_masks[class_index] = _haplotype_mask(
            class_haplotypes, haplotype_index
        )
        for product_id, sequence_hash, multiplicity in item.signature:
            product_position = product_index[product_id]
            key_position = exact_key_index[ProductKey(product_id, sequence_hash)]
            class_product[class_index, product_position] = 1
            class_exact[class_index, key_position] = 1
            class_multiplicity[class_index, product_position] += int(multiplicity)
        for product_id in _required_product_ids_for_class(space, item):
            if product_id in product_index:
                class_required[class_index, product_index[product_id]] = 1

    # key_index=-1 means no compiled-sequence hash was emitted. key_index=-2
    # means a hash was emitted but it is not predicted by any candidate class.
    # These two states have different likelihood semantics in the reference model.
    grouped: dict[tuple[int, int, int], list[float]] = defaultdict(
        lambda: [0.0, 0.0]
    )
    for item, weight in weighted_evidence:
        product_position = product_index[item.assigned_product_id]
        if item.best_sequence_sha256:
            key_position = exact_key_index.get(
                ProductKey(item.assigned_product_id, item.best_sequence_sha256),
                -2,
            )
        else:
            key_position = -1
        candidate_mask = int(
            _haplotype_mask(item.candidate_haplotypes, haplotype_index)
        )
        key = (product_position, key_position, candidate_mask)
        edit_penalty = (
            min(1.0, max(0.0, float(item.edit_rate)))
            if item.edit_rate is not None
            else 0.0
        )
        grouped[key][0] += float(weight)
        # The reference expression is weight * (log(probability) - edit_rate).
        grouped[key][1] += float(weight) * edit_penalty

    group_items = sorted(grouped.items())
    group_product = np.asarray(
        [key[0] for key, _values in group_items], dtype=np.int64
    )
    group_key = np.asarray(
        [key[1] for key, _values in group_items], dtype=np.int64
    )
    group_candidate_masks = np.asarray(
        [key[2] for key, _values in group_items], dtype=np.uint64
    )
    group_weight = np.asarray(
        [values[0] for _key, values in group_items], dtype=np.float64
    )
    group_edit_penalty = np.asarray(
        [values[1] for _key, values in group_items], dtype=np.float64
    )
    observed_counts = np.asarray(
        [float(observed.get(product_id, 0.0)) for product_id in product_ids],
        dtype=np.float64,
    )
    product_efficiencies = np.asarray(
        [float(config.product_efficiencies.get(product_id, 1.0)) for product_id in product_ids],
        dtype=np.float64,
    )
    priors = _class_priors(classes, config)
    log_priors = np.asarray(
        [log(priors[item.class_id]) for item in classes], dtype=np.float64
    )

    return (
        class_exact,
        class_product,
        class_required,
        class_multiplicity,
        class_haplotype_masks,
        group_product,
        group_key,
        group_candidate_masks,
        group_weight,
        group_edit_penalty,
        observed_counts,
        product_efficiencies,
        log_priors,
        useful_effective_reads,
    )


def score_genotypes_auto(
    space: CompiledGenotypeSpace,
    evidence: Sequence[MoleculeEvidence],
    config: GenotypeConfig | None = None,
    *,
    backend: BackendName = "auto",
) -> tuple[list[GenotypeClassScore], float, str]:
    """Score genotype classes with the compiled backend when it is safe.

    ``backend='python'`` is the semantic reference. ``backend='cython'`` requires
    the extension and raises rather than silently falling back. ``backend='auto'``
    uses Cython when available and the candidate catalogue fits the current dense
    uint64 haplotype-mask representation.
    """

    config = config or GenotypeConfig()
    config.validate()
    if backend not in {"auto", "python", "cython"}:
        raise GenotypeModelError(
            f"unknown genotype backend {backend!r}; expected auto, python or cython"
        )
    if backend == "python":
        scores, effective = score_genotypes(space, evidence, config)
        return scores, effective, "python_reference"
    if len(space.haplotype_ids) > 64:
        return _fallback_or_raise(
            backend=backend,
            reason="candidate catalogue contains more than 64 haplotypes",
            space=space,
            evidence=evidence,
            config=config,
        )
    try:
        from ._genotype_fast import score_dense_classes
    except ImportError:
        return _fallback_or_raise(
            backend=backend,
            reason="extension was not built",
            space=space,
            evidence=evidence,
            config=config,
        )

    classes = observable_classes(space)
    if not classes:
        raise GenotypeModelError("compiled assay yields no observable genotype classes")
    prepared = _prepare_dense_problem(
        space=space,
        classes=classes,
        evidence=evidence,
        config=config,
    )
    (
        class_exact,
        class_product,
        class_required,
        class_multiplicity,
        class_haplotype_masks,
        group_product,
        group_key,
        group_candidate_masks,
        group_weight,
        group_edit_penalty,
        observed_counts,
        product_efficiencies,
        log_priors,
        useful_effective_reads,
    ) = prepared
    read_ll, count_ll, dropout_ll, log_scores = score_dense_classes(
        class_exact,
        class_product,
        class_required,
        class_multiplicity,
        class_haplotype_masks,
        group_product,
        group_key,
        group_candidate_masks,
        group_weight,
        group_edit_penalty,
        observed_counts,
        product_efficiencies,
        log_priors,
        float(config.artifact_probability),
        float(config.product_background_mass),
        float(config.count_concentration),
        float(config.minimum_dirichlet_alpha),
        float(config.dropout_probability),
    )
    normalizer = _logsumexp([float(value) for value in log_scores])
    output = [
        GenotypeClassScore(
            class_id=item.class_id,
            genotype_pairs=item.genotype_pairs,
            signature=item.signature,
            log_read_likelihood=float(read_ll[index]),
            log_count_likelihood=float(count_ll[index]),
            log_dropout_likelihood=float(dropout_ll[index]),
            log_prior=float(log_priors[index]),
            log_score=float(log_scores[index]),
            posterior=exp(float(log_scores[index]) - normalizer),
        )
        for index, item in enumerate(classes)
    ]
    output.sort(key=lambda item: (-item.posterior, item.class_id))
    return output, useful_effective_reads, "cython_dense_grouped"


__all__ = [
    "BackendName",
    "cython_backend_available",
    "score_genotypes_auto",
]
