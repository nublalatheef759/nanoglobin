# NanoGlobin

NanoGlobin is an assay-aware research system for haemoglobinopathy analysis from
Oxford Nanopore reads and retrospective molecular-testing records. It targets the
`HBA1`/`HBA2` and beta-globin loci and treats the biological result as a pair of
chromosome haplotypes rather than a bag of unrelated SNV, indel and CNV calls.

The repository now contains four connected layers:

```text
retrospective reports
    -> provenance-preserving cohort/genotype ledger

candidate chromosome-haplotype sequences + declared primers
    -> in-silico PCR
    -> expected products + assay coverage + intrinsic ambiguity

FASTQ
    -> primer-aware molecule admission
    -> product, sequence and artifact evidence
    -> posterior over chromosome-haplotype pairs
    -> resolved / assay-equivalent / posterior-ambiguous / no-call

FASTQ or BAM
    -> Clair3 / bcftools csq / CIGAR / Sniffles / cuteSV / coverage evidence
    -> truth-scoped benchmarks and research reports
```

The generic callers remain useful evidence producers and baselines. The amplicon
model supplies semantics they do not own: primer observability, physical product
identity, PCR abundance, copy-marker ambiguity and assay-induced no-call states.

## Current programme

The useful boundary is **implementation versus validation**, not MSc versus
post-MSc.

Current method development includes:

- sequence-resolved chromosome-haplotype candidates;
- a generated HBA family-coordinate and copy-marker map;
- declared assay profiles and exact in-silico products;
- primer-aware FASTQ molecule classification;
- product and compiled-sequence compatibility;
- an initial posterior over pairs of chromosome haplotypes;
- configurable product efficiency, dropout and artifact terms;
- assay-equivalent ambiguity and explicit no-call states; and
- a Quarto/Typst thesis workspace that consumes pipeline results.

Analytical sensitivity, specificity, posterior calibration, lot/laboratory
robustness and clinical deployment claims remain gated by certified or
orthogonally characterised assay-matched material. Exact AmplideX support also
requires legally accessible assay documentation or reads.

See:

- [`docs/MSC_VERTICAL_SLICE.md`](docs/MSC_VERTICAL_SLICE.md) — active programme and validation gates;
- [`docs/HBA_FAMILY_MAP.md`](docs/HBA_FAMILY_MAP.md) — family coordinate and marker semantics;
- [`docs/GENOTYPE_MODEL.md`](docs/GENOTYPE_MODEL.md) — candidate posterior;
- [`docs/MOLECULE_ADMISSION.md`](docs/MOLECULE_ADMISSION.md) — primer-aware read evidence; and
- [`thesis/README.md`](thesis/README.md) — Quarto HTML and Typst PDF workflow.

## Core assay-aware workflow

### 1. Compile the assay

A profile declares primer sequences, product rules, length ranges, pools and
required products. Candidate chromosome haplotypes are supplied as FASTA
sequences. Each chromosome haplotype may encode one, two, three or another number
and order of alpha-like copies, hybrid genes and linked small variants.

```bash
python scripts/compile_assay.py \
  --profile assays/declared_assay.yml \
  --haplotypes catalogues/globin_haplotypes.fa \
  --compiled-json results/assay/compiled_products.json \
  --coverage-tsv results/assay/coverage.tsv \
  --indistinguishability-tsv results/assay/indistinguishability.tsv \
  --summary-json results/assay/summary.json
```

Outputs include exact product sequences and hashes, primer mismatches, product
multiplicity, assay coverage gaps and named genotype pairs that compile to the
same observable products.

### 2. Build the HBA family-coordinate map

The phrase “full HBA1/HBA2 population equivalence map” has been retired. The
implemented object is a generated **HBA family-coordinate and copy-marker map**.
It aligns HBA1, HBA2, hybrid and structural-copy sequences to one declared anchor,
retains substitutions/insertions/deletions, and reports copy-informative and
presence/absence markers.

```bash
python scripts/build_hba_family_map.py \
  --sequences catalogues/hba_copy_sequences.fa \
  --anchor-id HBA2_anchor \
  --positions-tsv results/family_map/positions.tsv \
  --markers-tsv results/family_map/markers.tsv \
  --map-json results/family_map/map.json \
  --summary-json results/family_map/summary.json
```

Adding HPRC or locally resolved haplotypes means appending sequences and
regenerating the map. It is not a manually curated all-population prerequisite.
Exact compiled product sequences remain the primary likelihood objects; the map
supports copy-marker interpretation, variant projection and visualization.

### 3. Admit and classify observed molecules

```bash
python scripts/admit_amplicon_reads.py \
  --sample SAMPLE \
  --fastq fastq/SAMPLE.fastq \
  --profile assays/declared_assay.yml \
  --compiled-json results/assay/compiled_products.json \
  --molecules-tsv results/assay/molecules/SAMPLE.molecules.tsv \
  --product-counts-tsv results/assay/molecules/SAMPLE.product_counts.tsv \
  --summary-json results/assay/molecules/SAMPLE.summary.json
```

Molecules are oriented from terminal primers and classified as complete,
one-ended, unexpected-length, primer-dimer, chimera-candidate or off-target.
Complete products are compared with compiled sequences while retaining equivalent
haplotypes, low-margin assignments and poor-fit residuals.

### 4. Rank chromosome-haplotype pairs

```bash
python scripts/genotype_amplicons.py \
  --compiled-json results/assay/compiled_products.json \
  --molecules-tsv results/assay/molecules/SAMPLE.molecules.tsv \
  --model-config examples/genotype_inference.example.yml \
  --posteriors-tsv results/assay/genotypes/SAMPLE.posteriors.tsv \
  --call-json results/assay/genotypes/SAMPLE.call.json \
  --summary-json results/assay/genotypes/SAMPLE.summary.json
```

The initial transparent model combines molecule/product compatibility, an
effective evidence cap for PCR non-independence, a Dirichlet-multinomial product
composition term, configurable product efficiencies, required-product dropout and
optional haplotype priors.

Named genotype pairs with the same count-aware product signature are collapsed
before the default prior is assigned. The output states are:

- `RESOLVED_RESEARCH_CALL`;
- `AMBIGUOUS_ASSAY_EQUIVALENT`;
- `AMBIGUOUS_POSTERIOR`; and
- `NO_CALL_INSUFFICIENT_EVIDENCE`.

These probabilities are conditional on the candidate catalogue, assay profile and
current parameters. They are explicitly labelled research/uncalibrated until
truth-matched material supports calibration.

## Legacy and baseline workflow

The established Snakemake path remains available:

```text
FASTQ -> minimap2 -> sorted/indexed BAM
      -> Clair3 phased SNV/indel candidates
      -> bcftools csq haplotype-aware consequence evidence
      -> CIGAR deletion evidence
      -> Sniffles and cuteSV structural candidates
      -> relative product/coverage profiles
      -> comprehensive and per-sample research reports
```

```bash
# Snakemake >= 8 and conda >= 24.7.1
snakemake --use-conda --cores 8 --rerun-triggers mtime
```

After changing a rule, put an explicit target after `--` when using
`--rerun-triggers mtime`:

```bash
snakemake --use-conda --cores 8 --rerun-triggers mtime -- <target>
```

NanoGlobin pins the existing workflow tools in `envs/`, including minimap2,
samtools, bcftools, Sniffles, cuteSV and Clair3. `config.yml` must point to the
Clair3 model matching the chemistry/basecaller.

## Evaluation evidence

Every result must be interpreted by assay, truth provenance, locus and variant
class. WGS comparisons do not validate multiplex-PCR behaviour, and agreement
with another caller is not independent truth.

### HPRC assembly-derived WGS comparison

Pipeline SNV/indel calls were compared with HPRC diploid-assembly callsets using
hap.py across five genomes:

| region | type | recall | precision |
|---|---|---:|---:|
| HBA/HBB core | SNV | 96.6% | 99.3% |
| HBA/HBB core | INDEL | 63.6% | 73.7% |
| beta-cluster gene bodies | SNV | 96.9% | 99.4% |
| beta-cluster gene bodies | INDEL | 63.3% | 91.2% |

Indels are the main limitation in this comparison. Assembly construction,
confidence scope and the exact included bases remain part of the result.

### GIAB HG002 regional comparison

Within released globin high-confidence regions and represented variant classes,
the current comparison reports 23/23 true positives, 0 false negatives and 0
false positives. This is a small regional small-variant result, not structural or
amplicon validation.

### HBA structural-event comparator set

Selected ONT samples were screened or labelled with the specialised DRAGEN HBA
caller and inspected with the CIGAR evidence path:

| sample | DRAGEN comparator | observed CIGAR evidence | zygosity estimate |
|---|---|---|---|
| NA21106 | -alpha4.2/alphaalpha | 4257 bp class | het (0.58) |
| HG00642 | -alpha3.7/alphaalpha | 3804 bp class | het |
| HG03136 | -alpha3.7/-alpha3.7 | 3804 bp class | hom (1.00) |
| HG00735 | alphaalphaalpha3.7/alphaalpha | no deletion event | excluded as deletion |
| HG03862 | --/alphaalpha | 10553 bp candidate | insufficient spanning evidence for zygosity |
| four comparator normals | alphaalpha/alphaalpha | none | — |

DRAGEN is a comparator/evidence producer, not independent ground truth. This
small comparator-selected set cannot estimate sensitivity or specificity.

Synthetic fixtures test known code paths and deliberate ambiguity. They are not
biological validation.

## CIGAR and coverage evidence

`scripts/detect_deletions_cigar.py` records large CIGAR `D` operations from
molecules bridging an event, compares broad size classes with `cnvs.csv`, and
reports breakpoint-spanning read fractions as zygosity evidence. The few-base
`-alpha3.7` subtype differences are not resolved by size alone.

Coverage profiles in amplicon mode are relative PCR-product abundance, not genomic
copy depth. They can flag unusual gain/loss shapes, but cannot alone establish
copy order, chromosome assignment, zygosity or a named triplication. Coverage-only
gains therefore remain unresolved structural candidates.

## Naming and annotation

The IthaGenes-primary catalogue (`databases/naming_layer.csv` converted to
`databases/variants.csv`) retains IthaGenes functionality, ClinVar assertions,
HbVar names and provenance. Unresolved coordinate derivations retain explicit
status rather than guessed positions.

IthaCNVs reciprocal-overlap and size matching provide candidate structural labels,
not exact chromosome haplotypes. Annotation is catalogue-guided when VEP returns
multiple transcripts, and `bcftools csq` supplies linked consequence evidence.

The REST annotation path remains convenient but mutable; a local versioned
annotation path is preferable for production-grade reproducibility.

## Retrospective cohort connection

The historical and ONT report tables contain selectively reported clinically
relevant findings, not complete detected callsets. Blank HBA/HBB fields are
therefore not converted to wild type. `normalise_cohort_reports.py` emits analysis
episodes, reported findings and phenotype measurements while retaining raw report
wording and explicit `reported`, `not_reported` and `not_tested` states.

The cohort work supplies the locally relevant haplotypes, HBA/HBB co-inheritance
states, phenotype questions, nomenclature failures and possible event-specific
truth sources. The assay model supplies an executable account of which of those
states a primer design can observe.

## Thesis workspace

`thesis/` is a Quarto book with shared sources for reviewable HTML and
submission-oriented Typst PDF. It includes chapter scaffolds for:

1. clinical and molecular background;
2. cohort design, ETL and data contracts;
3. reported genotype landscape and database representation;
4. HBA modification of HBB-associated phenotypes;
5. assay compiler, family coordinates and genotype model;
6. truth-scoped evaluation; and
7. discussion and next experiments.

Reusable R functions generate cohort-spectrum, assay-coverage, molecule-QC and
posterior figures. Chapters prefer approved pipeline results and fall back to
clearly labelled synthetic examples so the document remains renderable without
committing patient data.

```bash
cd thesis
quarto preview
quarto render --to html
quarto render --to typst
```

CI renders both formats and uploads the book artifact.

## Known limitations

- No currently accessible assay-matched truth cohort establishes analytical
  sensitivity or specificity for the public amplicon path.
- Product-efficiency and dropout defaults are explicit research parameters, not
  calibrated clinical constants.
- The family-coordinate implementation uses deterministic pairwise projection to
  one anchor; highly repetitive equally optimal alignments require unique-flank,
  graph or explicit candidate-haplotype adjudication.
- Automatic consensus/assembly of coherent unexplained products is not yet part of
  the merged workflow.
- Exact AmplideX integration requires permitted assay details or reads.
- Clair3 and generic SV caller outputs remain candidate evidence rather than the
  final assay-aware genotype.
- The current clinical-facing summaries still require consolidation around the
  chromosome-haplotype posterior and explicit no-call object.

## Repository layout

```text
Snakefile
config.yml
envs/

nanoglobin/
  assay.py                 assay profile and in-silico PCR
  family_map.py            HBA family coordinates and copy markers
  fastq.py                 streaming FASTQ/FASTQ.GZ input
  primer_matching.py       terminal primer recognition
  product_evidence.py      compiled product sequence evidence
  molecules.py             molecule admission orchestration
  molecule_reporting.py    molecule/product QC summaries
  genotype.py              chromosome-haplotype posterior
  cohort.py                selective-report cohort normalisation

scripts/
  compile_assay.py
  build_hba_family_map.py
  admit_amplicon_reads.py
  genotype_amplicons.py
  normalise_cohort_reports.py
  detect_deletions_cigar.py
  coverage_profile.py
  detect_cnv.py
  identify_sv.py
  vep_annotate.py
  comprehensive_report.py
  patient_summary.py
  sample_report.py

schemas/
examples/
databases/
tests/
docs/
thesis/
```
