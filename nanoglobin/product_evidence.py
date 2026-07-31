"""Compiled-assay loading and sequence evidence for amplicon molecules."""

from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from pathlib import Path
from typing import Iterable, Mapping, Sequence
import json
import math

from rapidfuzz.distance import Levenshtein

from .assay import (
    AssayProfile,
    AssayProfileError,
    CompiledProduct,
    normalise_dna,
)
from .molecule_types import (
    AdmissionConfig,
    CompiledAssay,
    MoleculeAdmissionError,
    SequenceCompatibility,
)

CompiledSequenceIndex = Mapping[
    str, Mapping[str, Sequence[CompiledProduct]]
]


def load_compiled_assay(
    path: str | Path,
    profile: AssayProfile | None = None,
) -> CompiledAssay:
    """Load and validate the compiler JSON used for read assignment."""

    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise MoleculeAdmissionError(
            f"cannot read compiled assay {path}: {exc}"
        ) from exc

    assay = payload.get("assay")
    raw_products = payload.get("compiled_products")
    if not isinstance(assay, Mapping) or not isinstance(raw_products, Sequence):
        raise MoleculeAdmissionError(
            "compiled assay JSON lacks assay or compiled_products"
        )

    assay_id = str(assay.get("assay_id", ""))
    version = str(assay.get("version", ""))
    if not assay_id or not version:
        raise MoleculeAdmissionError("compiled assay JSON lacks assay id/version")
    if profile and (assay_id != profile.assay_id or version != profile.version):
        raise MoleculeAdmissionError(
            f"compiled assay {assay_id}:{version} does not match profile "
            f"{profile.assay_id}:{profile.version}"
        )

    products: list[CompiledProduct] = []
    for index, raw_any in enumerate(raw_products):
        if not isinstance(raw_any, Mapping):
            raise MoleculeAdmissionError(
                f"compiled_products[{index}] is not an object"
            )
        raw = dict(raw_any)
        try:
            sequence = normalise_dna(
                str(raw["sequence"]),
                label=f"compiled product {index} sequence",
            )
            expected_hash = str(raw["sequence_sha256"])
            actual_hash = sha256(sequence.encode("ascii")).hexdigest()
            if expected_hash != actual_hash:
                raise MoleculeAdmissionError(
                    f"compiled product {index} sequence hash does not match sequence"
                )
            product = CompiledProduct(
                assay_id=str(raw["assay_id"]),
                assay_version=str(raw["assay_version"]),
                product_id=str(raw["product_id"]),
                pool=str(raw.get("pool", "default")),
                roles=tuple(str(value) for value in raw.get("roles", [])),
                required=bool(raw.get("required", False)),
                haplotype_id=str(raw["haplotype_id"]),
                start0=int(raw["start0"]),
                end0=int(raw["end0"]),
                length_bp=int(raw["length_bp"]),
                forward_mismatches=int(raw["forward_mismatches"]),
                reverse_mismatches=int(raw["reverse_mismatches"]),
                sequence=sequence,
                sequence_sha256=expected_hash,
            )
        except KeyError as exc:
            raise MoleculeAdmissionError(
                f"compiled product {index} lacks required field {exc.args[0]}"
            ) from exc
        except (TypeError, ValueError, AssayProfileError) as exc:
            raise MoleculeAdmissionError(
                f"invalid compiled product {index}: {exc}"
            ) from exc

        if product.assay_id != assay_id or product.assay_version != version:
            raise MoleculeAdmissionError(
                f"compiled product {index} does not share the declared assay id/version"
            )
        if product.length_bp != len(product.sequence):
            raise MoleculeAdmissionError(
                f"compiled product {index} length_bp does not match sequence length"
            )
        products.append(product)

    if not products:
        raise MoleculeAdmissionError("compiled assay contains no products")
    return CompiledAssay(
        assay_id=assay_id,
        assay_version=version,
        products=tuple(products),
    )


def build_compiled_sequence_index(
    products: Iterable[CompiledProduct],
) -> dict[str, dict[str, tuple[CompiledProduct, ...]]]:
    """Index equivalent compiled molecules by product identifier and sequence hash."""

    mutable: dict[str, dict[str, list[CompiledProduct]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for product in products:
        mutable[product.product_id][product.sequence_sha256].append(product)
    return {
        product_id: {
            sequence_hash: tuple(sorted(records, key=lambda item: item.haplotype_id))
            for sequence_hash, records in by_hash.items()
        }
        for product_id, by_hash in mutable.items()
    }


def score_compiled_sequences(
    observed: str,
    *,
    product_id: str,
    left_missing_bp: int,
    right_missing_bp: int,
    products: Iterable[CompiledProduct],
    config: AdmissionConfig,
    compiled_index: CompiledSequenceIndex | None = None,
) -> list[SequenceCompatibility]:
    """Score an admitted molecule against every compiled sequence for one product.

    The score is deliberately a simple edit-distance likelihood.  It is an
    inspectable first evidence model, not a calibrated ONT pair-HMM and not a
    diploid genotype posterior.
    """

    if compiled_index is None:
        by_hash = build_compiled_sequence_index(products).get(product_id, {})
    else:
        by_hash = compiled_index.get(product_id, {})

    output: list[SequenceCompatibility] = []
    error = config.assumed_sequence_error_rate
    for sequence_hash, records in by_hash.items():
        target_full = records[0].sequence
        target_end = (
            len(target_full) - right_missing_bp
            if right_missing_bp
            else len(target_full)
        )
        if left_missing_bp >= target_end:
            continue
        target = target_full[left_missing_bp:target_end]
        denominator = max(len(observed), len(target), 1)

        # RapidFuzz executes this dynamic programme in native code.  Use the
        # exact distance so a poor fit remains distinguishable from an absent
        # compiled sequence rather than being collapsed into no_compiled_sequence.
        distance = int(Levenshtein.distance(observed, target))
        edit_rate = distance / denominator
        log_likelihood = (
            (denominator - distance) * math.log(1 - error)
            + distance * math.log(error / 3)
        )
        output.append(
            SequenceCompatibility(
                product_id=product_id,
                sequence_sha256=sequence_hash,
                haplotype_ids=tuple(
                    sorted({record.haplotype_id for record in records})
                ),
                edit_distance=distance,
                edit_rate=edit_rate,
                log_likelihood=log_likelihood,
                target_length_bp=len(target),
                observed_length_bp=len(observed),
            )
        )

    return sorted(
        output,
        key=lambda item: (
            item.edit_rate,
            -item.log_likelihood,
            item.sequence_sha256,
        ),
    )
