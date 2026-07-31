"""Reusable NanoGlobin kernels."""

from .assay import (
    AssayProfile,
    AssayProfileError,
    CompiledProduct,
    Haplotype,
    assay_profile_from_mapping,
    compile_catalogue,
    coverage_rows,
    indistinguishability_classes,
    load_assay_profile,
    read_fasta,
)

__all__ = [
    "AssayProfile",
    "AssayProfileError",
    "CompiledProduct",
    "Haplotype",
    "assay_profile_from_mapping",
    "compile_catalogue",
    "coverage_rows",
    "indistinguishability_classes",
    "load_assay_profile",
    "read_fasta",
]

from .cohort import (
    CohortContract,
    CohortContractError,
    cohort_contract_from_mapping,
    load_cohort_contract,
    normalise_rows,
)

__all__ += [
    "CohortContract",
    "CohortContractError",
    "cohort_contract_from_mapping",
    "load_cohort_contract",
    "normalise_rows",
]
