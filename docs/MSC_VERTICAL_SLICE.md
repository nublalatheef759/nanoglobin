# MSc vertical slice and post-MSc programme

The repository contains a larger design programme because the correct biological
object is a pair of globin chromosome haplotypes and because an amplicon assay
observes primer-defined products rather than generic genomic depth. That does not
mean the entire programme must appear as completed MSc software.

The distinction below is based on scientific completeness, not on using the
calendar as an excuse for a generic pipeline.

The concept-review repository and unmerged agent branches are a design envelope,
not a list of features the student is expected to claim or reproduce. The
write-up should describe only merged, tested behaviour and should label any
remaining architecture as future work.

## Thesis-critical vertical slice

A coherent MSc software contribution is complete when one declared assay can be
followed through this chain:

```text
sequence-resolved candidate haplotypes
    -> in-silico PCR under a versioned assay profile
    -> assay coverage and intrinsic genotype ambiguity
    -> primer-aware admission of real ONT molecules
    -> product and artifact evidence
    -> sequence compatibility to compiled products
    -> explicit unresolved/no-call evidence
```

The thesis-critical implementation is therefore:

1. **Cohort contract and canonicalisation.** Preserve raw report wording and
   distinguish `reported`, `not_reported` and `not_tested`; do not manufacture
   wild-type calls from selective reports.
2. **One executable assay profile.** Use a primer design whose physical contract
   is available. A published profile can be used before private AmplideX data are
   accessible.
3. **Assay compiler.** Emit exact expected products, coverage gaps and genotype
   indistinguishability classes over a bounded globin haplotype catalogue.
4. **Molecule admission.** On real public amplicon FASTQ, report product
   assignment, read-length distributions, incomplete molecules, off-targets,
   primer dimers and chimera candidates.
5. **Sequence evidence.** Rank complete reads against compiled product sequences
   while preserving equivalent haplotypes, low-margin assignments and poor-fit
   residuals.
6. **Existing caller results as separate evidence.** Report the WGS small-variant
   benchmark and CIGAR/SV observations with their actual truth provenance; do not
   use them as assay-matched amplicon validation.
7. **Reproducible tests and visual QA.** Synthetic fixtures test code paths and
   known ambiguity; real unlabelled data test whether the assay model explains
   observed molecules.

This is already a real methods contribution. It is not merely a Snakemake wrapper:
the compiler and admission model encode primer observability, PCR molecule
identity and assay-induced ambiguity that generic callers do not represent.

## A defensible thesis result

The strongest result does not need to be a clinical sensitivity number. It can
be an executable assay audit:

- what fraction of public reads are explained by declared or inferred products;
- which products dominate or fail;
- observed product-length and orientation distributions;
- incomplete, off-target, dimer and chimera fractions;
- which candidate globin haplotypes the assay cannot distinguish by design;
- which reads remain coherent but unexplained by the catalogue; and
- how generic caller outputs agree or conflict with molecule-level evidence.

A result of "the assay cannot distinguish these states" is informative when it
is proven by compilation rather than hidden by a forced genotype.

## Optional extension if it becomes stable

A bounded diploid candidate ranker over known haplotypes is a reasonable stretch
result only after molecule admission is stable. It should expose sequence-only
and prior-free scores and must return ambiguity when product evidence is
insufficient.

It should not delay the cohort, assay-QC and write-up work merely to produce a
nominal genotype column.

## Post-MSc research programme

The following are extensions rather than prerequisites for a coherent thesis:

- a full HBA1/HBA2 sequence-equivalence map over population haplotypes;
- a variable-copy chromosome-haplotype posterior covering arbitrary HBA gains,
  hybrids and novel structures;
- hierarchical product-efficiency, dropout, lot and laboratory calibration;
- a calibrated ONT pair-HMM or neural read likelihood;
- automatic novel-product consensus and structural assembly;
- a legally supported AmplideX profile and raw-data integration;
- certified reference-material and blinded clinical validation;
- cross-assay and cross-laboratory portability; and
- clinical reporting or deployment claims.

Ground truth remains the final judge for those performance claims. The MSc slice
builds and falsifies the physical assay model so that eventual validation tests a
real method rather than a collection of correlated callers.

## Freeze criteria for the write-up

The implementation can be frozen for the thesis when:

- one profile and haplotype catalogue compile reproducibly;
- every input molecule is assigned to a product, artifact or unresolved state;
- assay-indistinguishable haplotypes remain indistinguishable;
- no blank report field becomes wild type;
- no coverage-only gain becomes a definitive named allele;
- public amplicon QC outputs are generated from exact declared samples;
- synthetic tests cover orientation, primer error, dropout, dimer, chimera,
  off-target and equivalent-haplotype cases; and
- all accuracy statements name their truth source, variant class and scope.

Everything beyond this gate can improve the project, but it is not needed to make
the software chapter scientifically coherent.
