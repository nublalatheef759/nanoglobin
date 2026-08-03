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

from .molecules import (
    AdmissionAccumulator,
    AdmissionConfig,
    AdmissionEngine,
    CompiledAssay,
    FastqRecord,
    MoleculeAdmissionError,
    MoleculeObservation,
    admit_molecule,
    load_compiled_assay,
    read_fastq,
)

__all__ += [
    "AdmissionAccumulator",
    "AdmissionConfig",
    "AdmissionEngine",
    "CompiledAssay",
    "FastqRecord",
    "MoleculeAdmissionError",
    "MoleculeObservation",
    "admit_molecule",
    "load_compiled_assay",
    "read_fastq",
]

from .family_map import (
    CoordinateObservation,
    FamilyMapError,
    MarkerObservation,
    SequenceRecord,
    align_to_anchor,
    build_coordinate_map,
    marker_rows,
)

__all__ += [
    "CoordinateObservation",
    "FamilyMapError",
    "MarkerObservation",
    "SequenceRecord",
    "align_to_anchor",
    "build_coordinate_map",
    "marker_rows",
]

from .genotype import (
    CompiledGenotypeSpace,
    GenotypeCall,
    GenotypeClassScore,
    GenotypeConfig,
    GenotypeModelError,
    MoleculeEvidence,
    analysis_payload,
    load_compiled_space,
    make_call,
    observable_classes,
    read_molecule_evidence,
    score_genotypes,
)

__all__ += [
    "CompiledGenotypeSpace",
    "GenotypeCall",
    "GenotypeClassScore",
    "GenotypeConfig",
    "GenotypeModelError",
    "MoleculeEvidence",
    "analysis_payload",
    "load_compiled_space",
    "make_call",
    "observable_classes",
    "read_molecule_evidence",
    "score_genotypes",
]

from .fast_genotype import (
    BackendName,
    cython_backend_available,
    score_genotypes_auto,
)

__all__ += [
    "BackendName",
    "cython_backend_available",
    "score_genotypes_auto",
]

from .calibration_types import (
    CalibrationError,
    CalibrationResult,
    CalibrationSettings,
    ComponentFit,
    ControlEvidence,
    ControlManifestRow,
)
from .calibration_fit import calibrate_model
from .calibration_io import (
    load_calibration_provenance,
    load_calibration_settings,
    read_control_manifest,
)

__all__ += [
    "CalibrationError",
    "CalibrationResult",
    "CalibrationSettings",
    "ComponentFit",
    "ControlEvidence",
    "ControlManifestRow",
    "calibrate_model",
    "load_calibration_provenance",
    "load_calibration_settings",
    "read_control_manifest",
]

from .truth_regions import (
    RegionAudit,
    TargetRegion,
    TruthRegionError,
    audit_regions,
    parse_target,
    read_bed,
)

__all__ += [
    "RegionAudit",
    "TargetRegion",
    "TruthRegionError",
    "audit_regions",
    "parse_target",
    "read_bed",
]
