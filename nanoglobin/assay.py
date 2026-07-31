"""Assay-profile validation and in-silico PCR compilation for NanoGlobin.

The module models the physical amplicon experiment before variant callers are
consulted.  A declared primer design is compiled against sequence-resolved
chromosome haplotypes.  The resulting product catalogue can be inspected for
coverage gaps and genotype pairs that are observationally indistinguishable.

The implementation intentionally uses only the Python standard library plus
PyYAML at the file-loading boundary.  Core functions accept ordinary mappings
and strings, which keeps them easy to test and reuse from Snakemake or another
workflow engine.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from hashlib import sha256
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence
import csv
import json
import re


IUPAC_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

IUPAC_COMPLEMENT: dict[str, str] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
}


class AssayProfileError(ValueError):
    """Raised when an assay profile is internally inconsistent."""


@dataclass(frozen=True)
class Primer:
    primer_id: str
    sequence: str
    max_mismatches: int = 0


@dataclass(frozen=True)
class ProductDefinition:
    product_id: str
    forward_primer: str
    reverse_primer: str
    min_length: int
    max_length: int
    pool: str = "default"
    roles: tuple[str, ...] = ()
    required: bool = False


@dataclass(frozen=True)
class AssayProfile:
    assay_id: str
    version: str
    platform: str
    primers: Mapping[str, Primer]
    products: tuple[ProductDefinition, ...]
    max_products_per_haplotype: int = 1_000


@dataclass(frozen=True)
class Haplotype:
    haplotype_id: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class PrimerHit:
    start: int
    end: int
    mismatches: int


@dataclass(frozen=True)
class CompiledProduct:
    assay_id: str
    assay_version: str
    product_id: str
    pool: str
    roles: tuple[str, ...]
    required: bool
    haplotype_id: str
    start0: int
    end0: int
    length_bp: int
    forward_mismatches: int
    reverse_mismatches: int
    sequence: str
    sequence_sha256: str

    def to_dict(self, *, include_sequence: bool = True) -> dict[str, Any]:
        value = asdict(self)
        value["roles"] = list(self.roles)
        if not include_sequence:
            value.pop("sequence")
        return value


@dataclass(frozen=True)
class IndistinguishabilityClass:
    class_id: str
    genotype_pairs: tuple[tuple[str, str], ...]
    product_signature: tuple[tuple[str, str, int], ...]


def normalise_dna(sequence: str, *, label: str = "sequence") -> str:
    """Return uppercase IUPAC DNA without whitespace, rejecting invalid bases."""

    value = re.sub(r"\s+", "", str(sequence)).upper().replace("U", "T")
    if not value:
        raise AssayProfileError(f"{label} is empty")
    invalid = sorted(set(value) - set(IUPAC_BASES))
    if invalid:
        raise AssayProfileError(
            f"{label} contains unsupported base(s): {', '.join(invalid)}"
        )
    return value


def reverse_complement(sequence: str) -> str:
    value = normalise_dna(sequence)
    return "".join(IUPAC_COMPLEMENT[base] for base in reversed(value))


def _bases_compatible(observed: str, pattern: str) -> bool:
    return bool(IUPAC_BASES[observed] & IUPAC_BASES[pattern])


def mismatch_count(observed: str, pattern: str) -> int:
    """Count IUPAC-aware substitutions between equal-length strings."""

    observed_n = normalise_dna(observed, label="observed sequence")
    pattern_n = normalise_dna(pattern, label="primer sequence")
    if len(observed_n) != len(pattern_n):
        raise ValueError("mismatch_count requires equal-length sequences")
    return sum(
        0 if _bases_compatible(obs, pat) else 1
        for obs, pat in zip(observed_n, pattern_n)
    )


def find_primer_hits(
    sequence: str,
    primer_sequence: str,
    *,
    max_mismatches: int,
) -> list[PrimerHit]:
    """Find all substitution-tolerant primer sites on a forward-oriented sequence.

    Indels at primer-binding sites are deliberately not hidden by this first
    compiler.  A profile can permit substitutions explicitly; an indel-bearing
    primer site remains non-amplifiable until a more expressive binding model is
    introduced and validated.
    """

    if max_mismatches < 0:
        raise ValueError("max_mismatches must be >= 0")
    target = normalise_dna(sequence)
    primer = normalise_dna(primer_sequence, label="primer sequence")
    width = len(primer)
    if width > len(target):
        return []

    hits: list[PrimerHit] = []
    for start in range(0, len(target) - width + 1):
        window = target[start : start + width]
        mismatches = sum(
            0 if _bases_compatible(observed, expected) else 1
            for observed, expected in zip(window, primer)
        )
        if mismatches <= max_mismatches:
            hits.append(PrimerHit(start=start, end=start + width, mismatches=mismatches))
    return hits


def _require_mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise AssayProfileError(f"{label} must be a mapping")
    return value


def _require_nonempty_string(value: Any, *, label: str) -> str:
    if value is None or not str(value).strip():
        raise AssayProfileError(f"{label} is required")
    return str(value).strip()


def assay_profile_from_mapping(raw: Mapping[str, Any]) -> AssayProfile:
    """Validate and construct an :class:`AssayProfile` from a mapping."""

    raw = _require_mapping(raw, label="assay profile")
    assay_id = _require_nonempty_string(raw.get("assay_id"), label="assay_id")
    version = _require_nonempty_string(raw.get("version", "1"), label="version")
    platform = _require_nonempty_string(raw.get("platform", "ont"), label="platform")

    primer_raw = _require_mapping(raw.get("primers"), label="primers")
    primers: dict[str, Primer] = {}
    for primer_id, spec_any in primer_raw.items():
        pid = _require_nonempty_string(primer_id, label="primer id")
        spec = _require_mapping(spec_any, label=f"primer {pid}")
        if pid in primers:
            raise AssayProfileError(f"duplicate primer id: {pid}")
        sequence = normalise_dna(spec.get("sequence", ""), label=f"primer {pid} sequence")
        max_mismatches = int(spec.get("max_mismatches", 0))
        if max_mismatches < 0:
            raise AssayProfileError(f"primer {pid} max_mismatches must be >= 0")
        if max_mismatches >= len(sequence):
            raise AssayProfileError(
                f"primer {pid} max_mismatches must be smaller than primer length"
            )
        primers[pid] = Primer(pid, sequence, max_mismatches)

    product_raw = raw.get("products")
    if not isinstance(product_raw, Sequence) or isinstance(product_raw, (str, bytes)):
        raise AssayProfileError("products must be a list")
    products: list[ProductDefinition] = []
    seen_products: set[str] = set()
    for index, item_any in enumerate(product_raw):
        item = _require_mapping(item_any, label=f"products[{index}]")
        product_id = _require_nonempty_string(
            item.get("id") or item.get("product_id"),
            label=f"products[{index}].id",
        )
        if product_id in seen_products:
            raise AssayProfileError(f"duplicate product id: {product_id}")
        seen_products.add(product_id)

        forward = _require_nonempty_string(
            item.get("forward_primer"), label=f"product {product_id} forward_primer"
        )
        reverse = _require_nonempty_string(
            item.get("reverse_primer"), label=f"product {product_id} reverse_primer"
        )
        missing = [pid for pid in (forward, reverse) if pid not in primers]
        if missing:
            raise AssayProfileError(
                f"product {product_id} references unknown primer(s): {', '.join(missing)}"
            )

        min_length = int(item.get("min_length", item.get("min_length_bp", 1)))
        max_length = int(item.get("max_length", item.get("max_length_bp", 100_000)))
        if min_length < 1 or max_length < min_length:
            raise AssayProfileError(
                f"product {product_id} has invalid length range {min_length}-{max_length}"
            )
        roles_any = item.get("roles", [])
        if isinstance(roles_any, str):
            roles = (roles_any,)
        elif isinstance(roles_any, Sequence):
            roles = tuple(str(role).strip() for role in roles_any if str(role).strip())
        else:
            raise AssayProfileError(f"product {product_id} roles must be a list or string")

        products.append(
            ProductDefinition(
                product_id=product_id,
                forward_primer=forward,
                reverse_primer=reverse,
                min_length=min_length,
                max_length=max_length,
                pool=str(item.get("pool", "default")),
                roles=roles,
                required=bool(item.get("required", False)),
            )
        )

    if not products:
        raise AssayProfileError("at least one product is required")
    max_products = int(raw.get("max_products_per_haplotype", 1_000))
    if max_products < 1:
        raise AssayProfileError("max_products_per_haplotype must be >= 1")

    return AssayProfile(
        assay_id=assay_id,
        version=version,
        platform=platform,
        primers=primers,
        products=tuple(products),
        max_products_per_haplotype=max_products,
    )


def load_assay_profile(path: str | Path) -> AssayProfile:
    """Load a JSON or YAML assay profile."""

    profile_path = Path(path)
    text = profile_path.read_text(encoding="utf-8")
    suffix = profile_path.suffix.lower()
    if suffix == ".json":
        raw = json.loads(text)
    else:
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover - dependency boundary
            raise RuntimeError("PyYAML is required to load YAML assay profiles") from exc
        raw = yaml.safe_load(text)
    return assay_profile_from_mapping(raw)


def read_fasta(path: str | Path) -> list[Haplotype]:
    """Read a small/medium FASTA catalogue without external dependencies."""

    records: list[Haplotype] = []
    current_id: str | None = None
    current_desc = ""
    chunks: list[str] = []

    def flush() -> None:
        nonlocal current_id, current_desc, chunks
        if current_id is None:
            return
        records.append(
            Haplotype(
                haplotype_id=current_id,
                description=current_desc,
                sequence=normalise_dna("".join(chunks), label=f"FASTA record {current_id}"),
            )
        )
        current_id = None
        current_desc = ""
        chunks = []

    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                if not header:
                    raise AssayProfileError(f"empty FASTA header at line {line_number}")
                parts = header.split(maxsplit=1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) == 2 else ""
            else:
                if current_id is None:
                    raise AssayProfileError(
                        f"FASTA sequence appears before the first header at line {line_number}"
                    )
                chunks.append(line)
    flush()

    if not records:
        raise AssayProfileError(f"FASTA catalogue {path} contains no records")
    id_counts = Counter(record.haplotype_id for record in records)
    duplicates = sorted(identifier for identifier, count in id_counts.items() if count > 1)
    if duplicates:
        raise AssayProfileError(f"duplicate FASTA identifiers: {', '.join(duplicates)}")
    return records


def compile_product_for_haplotype(
    profile: AssayProfile,
    product: ProductDefinition,
    haplotype: Haplotype,
) -> list[CompiledProduct]:
    """Compile one declared product against one chromosome haplotype."""

    sequence = normalise_dna(haplotype.sequence, label=haplotype.haplotype_id)
    forward = profile.primers[product.forward_primer]
    reverse = profile.primers[product.reverse_primer]
    reverse_binding = reverse_complement(reverse.sequence)

    forward_hits = find_primer_hits(
        sequence, forward.sequence, max_mismatches=forward.max_mismatches
    )
    reverse_hits = find_primer_hits(
        sequence, reverse_binding, max_mismatches=reverse.max_mismatches
    )

    compiled: list[CompiledProduct] = []
    for f_hit in forward_hits:
        for r_hit in reverse_hits:
            if r_hit.start < f_hit.end:
                continue
            length = r_hit.end - f_hit.start
            if not product.min_length <= length <= product.max_length:
                continue
            template_sequence = sequence[f_hit.start : r_hit.end]
            # PCR primers become part of the amplified molecule.  Template bases
            # at tolerated primer-site mismatches are therefore overwritten by
            # the declared oligonucleotide sequence after amplification.  Keeping
            # the raw template slice here would expose a base that the sequenced
            # amplicon does not physically contain and could falsely distinguish
            # haplotypes that differ only under a primer.
            interior_start = len(forward.sequence)
            interior_end = len(template_sequence) - len(reverse_binding)
            product_sequence = (
                forward.sequence
                + template_sequence[interior_start:interior_end]
                + reverse_binding
            )
            compiled.append(
                CompiledProduct(
                    assay_id=profile.assay_id,
                    assay_version=profile.version,
                    product_id=product.product_id,
                    pool=product.pool,
                    roles=product.roles,
                    required=product.required,
                    haplotype_id=haplotype.haplotype_id,
                    start0=f_hit.start,
                    end0=r_hit.end,
                    length_bp=length,
                    forward_mismatches=f_hit.mismatches,
                    reverse_mismatches=r_hit.mismatches,
                    sequence=product_sequence,
                    sequence_sha256=sha256(product_sequence.encode("ascii")).hexdigest(),
                )
            )
            if len(compiled) > profile.max_products_per_haplotype:
                raise AssayProfileError(
                    f"product {product.product_id} generated more than "
                    f"{profile.max_products_per_haplotype} products on "
                    f"{haplotype.haplotype_id}; tighten primer mismatch or length rules"
                )
    return compiled


def compile_catalogue(
    profile: AssayProfile,
    haplotypes: Iterable[Haplotype],
) -> list[CompiledProduct]:
    compiled: list[CompiledProduct] = []
    for haplotype in haplotypes:
        for product in profile.products:
            compiled.extend(compile_product_for_haplotype(profile, product, haplotype))
    return sorted(
        compiled,
        key=lambda item: (
            item.haplotype_id,
            item.product_id,
            item.start0,
            item.end0,
            item.forward_mismatches + item.reverse_mismatches,
        ),
    )


def coverage_rows(
    profile: AssayProfile,
    haplotypes: Iterable[Haplotype],
    compiled: Iterable[CompiledProduct],
) -> Iterator[dict[str, Any]]:
    by_key: dict[tuple[str, str], list[CompiledProduct]] = defaultdict(list)
    for product in compiled:
        by_key[(product.haplotype_id, product.product_id)].append(product)

    for haplotype in haplotypes:
        for definition in profile.products:
            products = by_key[(haplotype.haplotype_id, definition.product_id)]
            yield {
                "haplotype_id": haplotype.haplotype_id,
                "product_id": definition.product_id,
                "required": definition.required,
                "observable": bool(products),
                "n_products": len(products),
                "lengths_bp": ";".join(str(item.length_bp) for item in products),
                "forward_mismatches": ";".join(
                    str(item.forward_mismatches) for item in products
                ),
                "reverse_mismatches": ";".join(
                    str(item.reverse_mismatches) for item in products
                ),
            }


def _product_multiset_signature(
    products: Iterable[CompiledProduct],
) -> tuple[tuple[str, str, int], ...]:
    counts = Counter((item.product_id, item.sequence_sha256) for item in products)
    return tuple(
        sorted((product_id, sequence_hash, count) for (product_id, sequence_hash), count in counts.items())
    )


def indistinguishability_classes(
    haplotypes: Iterable[Haplotype],
    compiled: Iterable[CompiledProduct],
) -> list[IndistinguishabilityClass]:
    """Group diploid haplotype pairs with identical observable product multisets."""

    haplotypes_sorted = sorted(haplotypes, key=lambda item: item.haplotype_id)
    products_by_haplotype: dict[str, list[CompiledProduct]] = defaultdict(list)
    for product in compiled:
        products_by_haplotype[product.haplotype_id].append(product)

    groups: dict[
        tuple[tuple[str, str, int], ...], list[tuple[str, str]]
    ] = defaultdict(list)
    for left, right in combinations_with_replacement(haplotypes_sorted, 2):
        signature = _product_multiset_signature(
            [*products_by_haplotype[left.haplotype_id], *products_by_haplotype[right.haplotype_id]]
        )
        groups[signature].append((left.haplotype_id, right.haplotype_id))

    output: list[IndistinguishabilityClass] = []
    class_number = 1
    for signature, pairs in sorted(
        groups.items(), key=lambda item: (len(item[1]) * -1, item[1])
    ):
        if len(pairs) < 2:
            continue
        output.append(
            IndistinguishabilityClass(
                class_id=f"AIC{class_number:04d}",
                genotype_pairs=tuple(sorted(pairs)),
                product_signature=signature,
            )
        )
        class_number += 1
    return output


def write_compiled_json(
    path: str | Path,
    profile: AssayProfile,
    haplotypes: Iterable[Haplotype],
    compiled: Iterable[CompiledProduct],
) -> None:
    payload = {
        "schema_version": 1,
        "assay": {
            "assay_id": profile.assay_id,
            "version": profile.version,
            "platform": profile.platform,
        },
        "haplotypes": [
            {
                "haplotype_id": item.haplotype_id,
                "description": item.description,
                "length_bp": len(item.sequence),
                "sequence_sha256": sha256(item.sequence.encode("ascii")).hexdigest(),
            }
            for item in haplotypes
        ],
        "compiled_products": [item.to_dict(include_sequence=True) for item in compiled],
    }
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_coverage_tsv(path: str | Path, rows: Iterable[Mapping[str, Any]]) -> None:
    fieldnames = [
        "haplotype_id",
        "product_id",
        "required",
        "observable",
        "n_products",
        "lengths_bp",
        "forward_mismatches",
        "reverse_mismatches",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_indistinguishability_tsv(
    path: str | Path,
    classes: Iterable[IndistinguishabilityClass],
) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["class_id", "n_genotypes", "genotype_pairs", "product_signature"]
        )
        for item in classes:
            genotype_pairs = ";".join(f"{left}/{right}" for left, right in item.genotype_pairs)
            signature = ";".join(
                f"{product_id}:{sequence_hash}:{count}"
                for product_id, sequence_hash, count in item.product_signature
            )
            writer.writerow(
                [item.class_id, len(item.genotype_pairs), genotype_pairs, signature]
            )
