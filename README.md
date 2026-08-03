# NanoGlobin

NanoGlobin is an assay-aware research system for haemoglobinopathy analysis from
Oxford Nanopore reads and retrospective molecular-testing records. It targets the
`HBA1`/`HBA2` and beta-globin loci and treats the biological result as a pair of
chromosome haplotypes rather than a bag of unrelated SNV, indel and CNV calls.

The repository contains four connected layers:

```text
retrospective reports
    -> provenance-preserving cohort/genotype ledger

candidate chromosome-haplotype sequences + declared primers
    -> in-silico PCR
    -> expected products + assay coverage + intrinsic ambiguity

FASTQ
    -> primer-aware molecule admission
    -> product, sequence and artifact evidence
    -> calibrated posterior over chromosome-haplotype pairs
    -> resolved / assay-equivalent / posterior-ambiguous / no-call

FASTQ or BAM
    -> Clair3 / bcftools csq / CIGAR / Sniffles / cuteSV / coverage evidence
    -> truth-scoped benchmarks and research reports
```

The generic callers remain useful evidence producers and baselines. The amplicon
model supplies semantics they do not own: primer observability, physical product
identity, PCR abundance, copy-marker ambiguity, assay-induced no-call states and
explicit genotype indistinguishability.

## Current programme

The useful boundary is **implementation versus validation**, not MSc versus
post-MSc.

Current executable method development includes:

- sequence-resolved chromosome-haplotype candidates;
- a generated HBA family-coordinate and copy-marker map;
- declared assay profiles and exact in-silico products;
- primer-aware FASTQ molecule classification;
- product and compiled-sequence compatibility;
- a posterior over pairs of chromosome haplotypes;
- truth-scoped estimation of relative product efficiency, dropout and unsupported
  molecule mass;
- global and empirical-Bayes lot/run/operator calibration strata;
- assay-equivalent ambiguity and explicit no-call states;
- a grouped Cython genotype-scoring backend checked against the Python reference;
- exact benchmark-confidence-region auditing; and
- a Quarto/Typst thesis workspace that consumes pipeline results.

Analytical sensitivity, specificity, posterior calibration, lot/laboratory
robustness and clinical deployment claims remain gated by certified or
orthogonally characterised assay-matched material. Exact AmplideX support also
requires legally accessible assay documentation or reads.

See:

- [`docs/MSC_VERTICAL_SLICE.md`](docs/MSC_VERTICAL_SLICE.md) — active programme and validation gates;
- [`docs/HBA_FAMILY_MAP.md`](docs/HBA_FAMILY_MAP.md) — family-coordinate and marker semantics;
- [`docs/GENOTYPE_MODEL.md`](docs/GENOTYPE_MODEL.md) — candidate posterior;
- [`docs/MODEL_CALIBRATION.md`](docs/MODEL_CALIBRATION.md) — calibration, identifiability and Cython execution;
- [`docs/MOLECULE_ADMISSION.md`](docs/MOLECULE_ADMISSION.md) — primer-aware read evidence;
- [`validation/README.md`](validation/README.md) — validation-data acquisition and control design; and
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

The repository includes
[`assays/huang_2023_ont_long_pcr.yml`](assays/huang_2023_ont_long_pcr.yml), an
executable transcription of a published multiplex long-PCR design. It is a public
assay analogue and **not** an AmplideX profile.

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

Full-sequence edit distance already executes in RapidFuzz's native kernel. The
Cython work therefore targets the remaining class-by-evidence genotype scorer,
not code that is already native.

### 4. Calibrate assay-product terms

A tab-separated control manifest declares the known chromosome-haplotype pair,
truth source, assay scope and calibration stratum for every control. WGS and
caller-comparator controls are excluded by default because they do not estimate
multiplex-PCR efficiency or independent truth.

```bash
python scripts/calibrate_amplicon_model.py \
  --compiled-json results/assay/compiled_products.json \
  --control-manifest validation/controls.tsv \
  --base-model-config examples/genotype_inference.example.yml \
  --calibration-settings examples/calibration_settings.example.yml \
  --calibration-json results/assay/calibration/calibration.json \
  --calibrated-model-config results/assay/calibration/calibrated_model.yml \
  --product-parameters-tsv results/assay/calibration/product_parameters.tsv \
  --sample-product-fit-tsv results/assay/calibration/sample_product_fit.tsv
```

The fitter estimates:

- relative product efficiency within identifiable connected product groups;
- required-product dropout with explicit Beta priors;
- unsupported complete-molecule mass;
- optional run/lot/operator strata shrunk toward the global fit; and
- every excluded control and identifiability warning.

One isolated product has no independently identifiable relative efficiency, and
disconnected product groups have no common observed scale. Those states are
reported rather than inferred from total library size.

### 5. Rank chromosome-haplotype pairs

```bash
python scripts/genotype_amplicons.py \
  --compiled-json results/assay/compiled_products.json \
  --molecules-tsv results/assay/molecules/SAMPLE.molecules.tsv \
  --model-config results/assay/calibration/calibrated_model.yml \
  --calibration-json results/assay/calibration/calibration.json \
  --backend auto \
  --posteriors-tsv results/assay/genotypes/SAMPLE.posteriors.tsv \
  --call-json results/assay/genotypes/SAMPLE.call.json \
  --summary-json results/assay/genotypes/SAMPLE.summary.json
```

The transparent model combines molecule/product compatibility, an effective
evidence cap for PCR non-independence, a Dirichlet-multinomial product-composition
term, product efficiencies, required-product dropout and optional haplotype priors.

Named genotype pairs with the same count-aware product signature are collapsed
before the default prior is assigned. The output states are:

- `RESOLVED_RESEARCH_CALL`;
- `AMBIGUOUS_ASSAY_EQUIVALENT`;
- `AMBIGUOUS_POSTERIOR`; and
- `NO_CALL_INSUFFICIENT_EVIDENCE`.

The calibration document embeds the exact fitted model configuration. Inference
fails if a different config is supplied under the same calibration provenance.
The probabilities remain conditional on the candidate catalogue, assay profile
and declared model.

### 6. Accelerate and benchmark the posterior

Editable/standard package installs build the Cython extension through PEP 517.
For a source checkout:

```bash
python setup.py build_ext --inplace
```

`--backend python` is the semantic reference. `--backend cython` fails closed if
the extension is unavailable. `--backend auto` uses the grouped dense Cython
kernel when the current candidate catalogue is representable and otherwise
records the explicit Python fallback reason.

```bash
python scripts/benchmark_genotype_backend.py \
  --compiled-json results/assay/compiled_products.json \
  --molecules-tsv results/assay/molecules/SAMPLE.molecules.tsv \
  --model-config results/assay/calibration/calibrated_model.yml \
  --warmups 2 --repeats 10 \
  --output-json results/benchmarks/SAMPLE.genotype_backend.json
```

The benchmark requires numerical parity across all likelihood terms and posterior
probabilities before reporting every repeated timing, median, spread and speedup.
No speed claim is made from a single favourable run.

A standalone calibrated Snakemake layer is available:

```bash
snakemake -s Snakefile.calibration --use-conda --cores 4
```

### 7. Audit benchmark truth coverage

A WGS benchmark applies only where its confidence BED covers the exact locus and
variant class. Audit the intervals before running or interpreting `hap.py`:

```bash
python scripts/audit_truth_regions.py \
  --bed truth/HG002.confident.bed \
  --target HBA=chr16:170000-178000 \
  --target HBB=chr11:5225000-5310000 \
  --output-tsv results/truth/HG002.globin_coverage.tsv \
  --summary-json results/truth/HG002.globin_coverage.json
```

`--require-full-coverage` returns a distinct non-zero status when any requested
interval is incomplete.

## Validation data programme

[`validation/resources.tsv`](validation/resources.tsv) is the executable resource
registry. The acquisition order is:

1. use public multiplex-PCR reads for molecule mechanics, endpoint reconstruction,
   artifact distributions and runtime—not genotype accuracy without truth;
2. audit HPRC and GIAB at sequence and confidence-region level for candidate
   haplotypes and WGS caller behaviour;
3. request raw reads and sample-level truth from published truth-matched ONT
   haemoglobinopathy cohorts;
4. order and orthogonally reconfirm available Coriell/CDC HBA and HBB controls;
5. recover an exact AmplideX blinded set with the historical pre-ONT GAP-PCR,
   MLPA or Sanger result for each sample; and
6. fit on a training/control subset, freeze the assay profile, catalogue and model,
   then evaluate an independent held-out set.

Vendor Reporter and DRAGEN outputs remain comparators unless the genotype is
independently adjudicated. Validation must retain negatives, failures, ambiguous
results and no-calls, and must be stratified by event class.

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
amplicon validation. The exact confidence BED should be re-audited with
`audit_truth_regions.py` for every frozen benchmark release.

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
5. assay compiler, family coordinates, calibration and genotype model;
6. truth-scoped evaluation and performance; and
7. discussion and next experiments.

Reusable R functions generate cohort-spectrum, assay-coverage, molecule-QC,
calibration, posterior and backend-performance figures. Chapters prefer approved
pipeline results and fall back to clearly labelled synthetic examples so the
document remains renderable without committing patient data.

```bash
cd thesis
quarto preview
quarto render --to html
quarto render --to typst
```

CI renders both formats and uploads the book artifact.

## Known limitations

- No currently accessible exact-assay truth cohort establishes analytical
  sensitivity or specificity for the AmplideX path.
- Control-fitted parameters remain conditional on the controls, event classes,
  runs and candidate catalogue used; they are not universal PCR constants.
- The current likelihood uses one global required-product dropout probability for
  inference while retaining product-specific dropout estimates as diagnostics.
- The current Cython dense backend uses a `uint64` haplotype mask and falls back to
  the Python reference for catalogues larger than 64 named candidates.
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
Snakefile.calibration
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
  genotype.py              Python reference chromosome-haplotype posterior
  fast_genotype.py         automatic grouped Python/Cython scoring
  _genotype_fast.pyx       typed dense scoring kernel
  calibration_*.py         truth-scoped calibration and provenance
  truth_regions.py         benchmark-confidence-region auditing
  cohort.py                selective-report cohort normalisation

scripts/
  compile_assay.py
  build_hba_family_map.py
  admit_amplicon_reads.py
  calibrate_amplicon_model.py
  genotype_amplicons.py
  benchmark_genotype_backend.py
  audit_truth_regions.py
  normalise_cohort_reports.py
  detect_deletions_cigar.py
  coverage_profile.py
  detect_cnv.py
  identify_sv.py
  vep_annotate.py
  comprehensive_report.py
  patient_summary.py
  sample_report.py

assays/
validation/
schemas/
examples/
databases/
tests/
docs/
thesis/
```
