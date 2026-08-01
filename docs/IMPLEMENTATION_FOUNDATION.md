# Cohort-to-caller implementation foundation

This changeset introduces two executable contracts that connect the retrospective
cohort work to the future amplicon genotyper.

## 1. Cohort report normalisation

`normalise_cohort_reports.py` converts a wide, selectively reported clinical
spreadsheet into three long-form tables:

- analysis episodes;
- reported locus findings; and
- phenotype measurements.

It preserves the raw report wording and never turns an empty HBA/HBB cell into a
reference genotype. Missing text is `not_reported`; an explicit sentinel is
`not_tested`. The generated IDs are provenance identifiers, not inferred patient
linkage. A real de-identified person identifier should be declared in the cohort
contract when available.

```bash
python scripts/normalise_cohort_reports.py \
  --input governed/historical.csv \
  --contract governed/historical.contract.yml \
  --episodes results/cohort/episodes.csv \
  --findings results/cohort/findings.csv \
  --phenotypes results/cohort/phenotypes.csv \
  --summary-json results/cohort/summary.json
```

The example contract under `examples/` contains column names only and no patient
data.

## 2. Assay compilation

`compile_assay.py` validates a versioned primer/product profile and runs in-silico
PCR against chromosome-haplotype FASTA sequences. It emits:

- exact compiled product sequences and hashes;
- an assay coverage table showing missing required products; and
- diploid genotype pairs with identical observable product multisets.

```bash
python scripts/compile_assay.py \
  --profile assays/huang_2023.yml \
  --haplotypes catalogues/globin_haplotypes.fa \
  --compiled-json results/assay/compiled_products.json \
  --coverage-tsv results/assay/coverage.tsv \
  --indistinguishability-tsv results/assay/indistinguishability.tsv \
  --summary-json results/assay/summary.json
```

The compiler does not call genotypes. It establishes what the assay can
physically amplify and which genotype states cannot be separated even with
perfect reads. Primer-site substitutions are supported through an explicit
mismatch allowance. Primer-site indels remain non-amplifiable in this first,
auditable model rather than being silently approximated.

## 3. Report provenance and structural-gain semantics

Clinical-facing reports now receive the configured sample list explicitly.
Filesystem globbing is removed from `comprehensive_report.py` and `detect_cnv.py`,
so stale or simulated samples cannot enter a report merely because their files
exist.

Coverage-defined HBA gains are reported as structural candidates. Reciprocal
interval overlap may identify a catalogue comparator, but it does not establish
copy order, junction sequence, or a named allele such as an anti-3.7
triplication. `sample_report.py` therefore ranks an unconfirmed gain as an
unresolved finding rather than a pathogenic call.

## 4. Molecule admission

The second executable layer is now implemented in
[`MOLECULE_ADMISSION.md`](MOLECULE_ADMISSION.md). It recognises terminal primers,
orients reads, preserves incomplete/dimer/chimera/off-target states and scores
complete molecules against compiled product sequences. It does not force a
diploid genotype.

The bounded MSc vertical slice and the distinction from the longer research
programme are defined in [`MSC_VERTICAL_SLICE.md`](MSC_VERTICAL_SLICE.md).

## 5. What remains

The next engine slices are a HBA sequence-equivalence map, product-abundance and
dropout models, and a posterior over chromosome-haplotype pairs. Analytical
sensitivity and specificity remain contingent on sealed, assay-matched ground
truth.
