# NanoGlobin

Variant calling for haemoglobinopathies from Oxford Nanopore long-read data —
α/β-thalassaemias, structural haemoglobin variants (Hb S, Hb C, Hb E), and
copy-number rearrangements. Targets HBA1/HBA2 (chr16) and HBB (chr11);
developed for thalassaemia (nanothal) and generalised to the globin loci
(nanoglobin).

```
assay profile + haplotype FASTA → compile_assay → expected products + intrinsic ambiguity
FASTQ → admit_amplicon_reads → molecule/product assignment + artifact/QC evidence
FASTQ → minimap2 → sort/index ─┬→ Clair3 (SNV/indel, phased) → filter → region filter → VEP → annotate
                               ├→ bcftools csq (haplotype-aware consequence)
                               ├→ CIGAR deletion evidence
                               ├→ Sniffles ┐
                               ├→ CuteSV   ┤→ identify_sv (IthaCNVs comparator)
                               └→ coverage_profile (relative product abundance)
                                            └→ comprehensive report → patient summary → research report
```

## Status and scope

NanoGlobin is a research pipeline with component-level WGS benchmarks,
comparator checks for selected HBA structures, and an assay-aware amplicon
evidence layer. These are complementary evidence sources, not three independent
truth sets and not assay-matched clinical validation.

The merged MSc vertical slice and the larger post-MSc programme are separated in
[`docs/MSC_VERTICAL_SLICE.md`](docs/MSC_VERTICAL_SLICE.md). Primer-aware read
admission is documented in
[`docs/MOLECULE_ADMISSION.md`](docs/MOLECULE_ADMISSION.md). Only merged and tested
behaviour should be described as implemented; the larger design is future work.

The pipeline runs anywhere via `--use-conda` (see Reproducibility). The legacy
BlueBEAR module path still works but is not required.

> **Running the pipeline after adding/editing a rule:** Snakemake treats a changed
> rule as a reason to rebuild everything upstream. Always pass
> `--rerun-triggers mtime` so only genuinely stale targets rebuild, and put the
> target *after* a `--` separator (the flag otherwise swallows it):
> `snakemake --use-conda --cores 8 --rerun-triggers mtime -- <target>`.

## Reproducibility

Every rule runs in its own pinned conda environment (`envs/*.yaml`), created
automatically by Snakemake. There is no module loading and no shell-switching —
the whole pipeline runs from a single command.

```bash
# one-time: snakemake (>=8) with conda >=24.7.1 in its environment

# one-time: obtain the Clair3 ONT model (not bundled with the conda package).
# NanoGlobin pins Clair3 1.0.4, which uses v1 (TensorFlow) models. Download
# r1041_e82_400bps_sup_v430 — or the model matching your basecaller — from ONT
# rerio (https://github.com/nanoporetech/rerio), place it under models/, and set
# config.yml `clair3_models:` to that path (e.g. models/r1041_e82_400bps_sup_v430).

# run everything
snakemake --use-conda --cores 8 --rerun-triggers mtime
```

`config.yml` `clair3_models:` must point at the extracted model directory. The
`r1041_e82_400bps_sup_v430` model matches R10.4.1 E8.2 (5kHz) Dorado v4.3.0 SUP
basecalling; data from a different chemistry/basecaller needs the corresponding
model from rerio. rerio Clair3 models are v1-compatible — do not upgrade to
Clair3 2.x without a matching v2 model.

Pinned tool versions (`envs/`): minimap2 2.24, samtools 1.16.1, bcftools 1.15.1,
sniffles 2.8.0, cuteSV 1.0.8, Clair3 1.0.4. The SV-caller environments list
conda-forge before bioconda so pysam/libdeflate resolve under strict channel
priority; the cuteSV rule sets `PYTHONNOUSERSITE=1` so a stray user-site cuteSV
install cannot shadow the pinned one.

## Evaluation evidence

The current evidence must be interpreted by assay, truth provenance, locus and
variant class. WGS comparisons do not validate multiplex-PCR behaviour, and a
comparison with another caller is not independent truth.

**1. SNV/indel — HPRC assembly-derived WGS comparison (5 genomes).**
Pipeline calls were compared with HPRC diploid-assembly callsets using hap.py
(jmcdani20/hap.py v0.3.12), pooled across five genomes (HG02071, HG02083,
HG02514, HG02074, HG02622):

| region | type | recall | precision |
|---|---|---|---|
| HBA/HBB core | SNV | 96.6% | 99.3% |
| HBA/HBB core | INDEL | 63.6% | 73.7% |
| β-cluster gene bodies | SNV | 96.9% | 99.4% |
| β-cluster gene bodies | INDEL | 63.3% | 91.2% |

SNV calling is strong across both the core loci and the extended β-cluster gene
bodies in this comparison. Indels are the primary limitation (~63% recall), with
false negatives clustering at a recurring low-complexity tandem-repeat region
(chr16:171,206–171,221). Across the full extended interval (including intergenic
and locus-control-region sequence), SNV precision falls to 94.1%, driven by
false positives in repetitive/regulatory intergenic and LCR sequence rather than
in the globin genes themselves. Assembly-derived callsets and confidence scope
must remain explicit when these numbers are reported.

**2. SNV/indel — GIAB HG002 within declared confident regions.**
Against NIST v4.2.1, restricted to globin high-confidence regions and represented
variant classes: **23/23 true positives, 0 FN, 0 FP = 100% / 100%.** This is a
small clean regional benchmark, not a locus-complete structural or amplicon
validation result.

**3. HBA structural-event comparator set (CIGAR method).**
Selected ONT samples were screened or labelled with the specialised DRAGEN HBA
caller and then inspected with the CIGAR method. DRAGEN is a comparator/evidence
producer here, not independent ground truth. The small selected set showed no
false-positive deletion calls among the included controls, but it cannot estimate
sensitivity or specificity.

| sample | DRAGEN comparator | detected | zygosity estimate |
|---|---|---|---|
| NA21106 | -α4.2/αα | 4257 bp = -α4.2 | het (0.58) |
| HG00642 | -α3.7/αα | 3804 bp = -α3.7 | het |
| HG03136 | -α3.7/-α3.7 | 3804 bp = -α3.7 | hom (1.00) |
| HG00735 | ααα3.7/αα | no deletion event | excluded as deletion |
| HG03862 | --/αα | 10553 bp candidate | insufficient spanning evidence for zygosity |
| normals ×4 | αα/αα | none | — |

Synthetic positives (`scripts/simulate/`) are deterministic development fixtures
for code paths and known-answer cases. They are not biological validation.

## α-globin deletion evidence (CIGAR method)

`scripts/detect_deletions_cigar.py`. In the current alignment representation,
ONT reads spanning a known α-deletion can carry a large CIGAR `D`. The method:

- records large deletion evidence from spanning-read CIGAR operations;
- compares the observed size with `cnvs.csv`, separating the broad `-α3.7` and
  `-α4.2` classes in the current examples;
- reports the deletion-read fraction at the breakpoint as a zygosity estimate;
  and
- does not convert absence of a deletion CIGAR into a triplication call.

The `-α3.7` type I/II/III subtypes differ by only a few bases and are not resolved
by breakpoint size. Class-level reporting is used.

This supplies direct molecule-level junction evidence when a read bridges the
event, but the present set is too small and comparator-selected to establish
analytical accuracy. Large events with insufficient spanning molecules remain
unresolved by this path.

**Coverage evidence** (`scripts/coverage_profile.py`, `detect_cnv.py`) is binned
and currently normalised to the configured HBB product/region. In an amplicon
assay, this is relative PCR-product abundance rather than genomic copy depth. It
can flag unusual loss or gain shapes for review, but cannot by itself establish
HBA copy number, zygosity, copy order or a named triplication. The report therefore
retains coverage-defined gains as unresolved structural candidates.

## Assay compilation and molecule admission

When `assay_profile` and `haplotype_catalogue` are configured, NanoGlobin can:

1. compile primer-defined products over sequence-resolved candidate haplotypes;
2. identify assay coverage gaps and intrinsic genotype indistinguishability;
3. recognise terminal primers in both FASTQ orientations;
4. classify complete, one-ended, dimer, unexpected-length, chimera and off-target
   molecules; and
5. compare complete molecules with compiled product sequences while retaining
   equivalent-haplotype and low-margin assignments.

```bash
snakemake --use-conda --cores 8 \
  results/assay/molecules/SAMPLE.summary.json
```

These outputs are assay-product and artifact evidence. Product counts are not
genomic allele counts, and this layer deliberately stops before a diploid
chromosome-haplotype posterior.

## Naming and classification layer

Variants are named and classified against an **IthaGenes-primary catalogue**
(`databases/naming_layer.csv` → converted to `databases/variants.csv`, 3,685
entries), keyed on HGVS (`GENE:c.notation`). This supersedes the earlier
ClinVar-derived catalogue (2,839 entries) — it is a strict superset (+846
variants) that uses expert thalassaemia curation as the classification spine:

- **IthaGenes (ITHANET)** — 2,081 HBA/HBB variants (2,028 Causative + 53 Neutral),
  the curation spine. `scripts/build_naming_layer.py`.
- **ClinVar** — 1,604 additional variants not in IthaGenes, plus classification
  cross-reference where both hold a variant.
- Source per variant: 1,069 IthaGenes-only, 1,012 IthaGenes+ClinVar,
  1,604 ClinVar-only.

Classification is **Functionality-first**: `sample_report.py:tier()` checks the
IthaGenes Functionality (`Causative`) before ClinVar, so IthaGenes causative
variants tier as pathogenic even where ClinVar is silent or VUS-heavy. No
ClinVar terms are fabricated for IthaGenes-only variants; a blank ClinVar column
remains visible.

**GRCh38 coordinates** were derived via the VariantValidator REST API and checked
against 123 known coordinates (45 general + 78 HBB promoter) before application
to unresolved catalogue entries. Coverage is reported as **99.8% (3,676/3,685)**.
Derivation paths include:

- direct GRCh38 VCF blocks for in-transcript variants;
- genomic re-submit suggestions for out-of-transcript substitutions;
- transcript-to-genomic queries for out-of-transcript indels;
- component splitting for compound `[a;b]` alleles; and
- explicit HGVS cleaning rules.

The nine unresolved entries retain a `coord_status` reason and are not force-filled
with guessed coordinates.

Structural events are compared with **IthaCNVs** (`databases/cnvs.csv`, 311 CNVs,
GRCh38.p13). Reciprocal overlap and size can provide a candidate catalogue class,
but do not uniquely establish a chromosome haplotype. Coverage-defined gains are
not promoted to named triplications without junction/copy-specific evidence.

Annotation is **catalogue-guided**: where VEP returns several transcripts, the one
whose HGVS matches a catalogue entry is preferred, resolving promoter/5′UTR
variants that canonical-only annotation leaves as `upstream_gene_variant`.
Haplotype-aware consequences are additionally produced by **bcftools csq**
(`rule csq_annotate`, Ensembl GFF3), which reports linked consequences rather than
only isolated per-site consequences.

## Co-inheritance flagging

Patients co-inheriting **causative** HBA and HBB findings are flagged for review.
In the parsed historical referral cohort, the fraction of reported HBB
heterozygotes below the 3.5% HbA2 cutoff rose across the current HBA report groups:
25.5% (no reported HBA variant), 36.5% (HBA heterozygous), and 43.3% (HBA
homozygous/compound heterozygous). These are retrospective referral-cohort
associations, not population estimates or causal validation. The analysis should
be regenerated from the governed cohort contract with explicit denominators,
missingness and covariates before thesis freeze.

## Per-sample research report

`results/sample_report.csv` gives one row per sample for high-volume review: a
descriptive Result flag, primary finding, HGVS and zygosity, counts by significance
tier, and the full ranked finding list. It is a research summary, not an automated
diagnosis.

SNV and structural-event evidence are unified into one view. Coverage-only gain
candidates remain explicitly unresolved. Benign and VUS findings sort below
causative findings but remain visible; conflicting classifications are surfaced
explicitly.

## Variant filtering and phasing

**Filtering** — in the planted Hb S fixture, `HBB:c.20A>T` at DP 991 and AF 0.465
was dropped by the previous `QUAL>20` rule because QUAL was 19.88, while low-depth
artifacts could retain high apparent AF. The current fixture-driven filter is
`FILTER=PASS` + `FORMAT/DP≥10` + `FORMAT/AF≥0.15`. These thresholds remain
candidate filters and require assay/chemistry-specific validation.

**Phasing** — Clair3 runs with `--enable_phasing`. In SRR37686273, six reported
heterozygous variants across 1.3 kb of HBB were emitted in one phased
configuration and can be inspected against individual molecules.

## Known issues

**Callers**
- CuteSV reports `./.` in current examples even with `--genotype`; its missing
  genotype is not used as authoritative zygosity.
- Clair3 requires the target BED plus downstream region filtering; the existing
  generic caller path is separate from primer-aware molecule admission.

**Annotation**
- Promoter/5′UTR HGVS is transcript-dependent; catalogue-guided selection is
  required. Variants not in the catalogue and sufficiently far upstream may
  remain `upstream_gene_variant`.
- `vep_annotate.py` uses the Ensembl REST region endpoint and can return an API
  error for some insertion representations. A local versioned annotation path is
  preferable for a production workflow.

**Clinical logic**
- Historical β⁰/β⁺ classes remain a curated data file
  (`beta_classification.csv`) and require provenance/version review.
- A catalogue match names an annotation candidate; it does not replace event
  sequence, phase, dosage evidence or clinical adjudication.

**Pipeline**
- WGS mode switches the reference, but no complete WGS cohort has been run
  end-to-end through every reporting layer.
- `patient_summary.py` consumes the explicitly constructed comprehensive report;
  report semantics still need consolidation around chromosome-haplotype and
  unresolved/no-call objects.

## Data sources

- **IthaGenes / IthaCNVs** (https://www.ithanet.eu) — variant curation, common
  names, functionality, and CNV breakpoints (GRCh38.p13). `build_naming_layer.py`,
  `refresh_cnvs.py`, `parse_ithacnv.py`.
- **ClinVar** — `variant_summary.txt.gz`, merged by `scripts/build_clinvar.py`.
- **HbVar** (https://globin.bx.psu.edu/hbvar) — common names, merged by
  `scripts/merge_hbvar.py`.
- **VariantValidator** (https://rest.variantvalidator.org) — GRCh38 coordinate
  derivation (MANE transcripts NM_000518.5 HBB, NM_000558.5 HBA1,
  NM_000517.6 HBA2).

## Layout

```
Snakefile                        pipeline definition and optional assay-evidence rules
config.yml                       paths, candidate thresholds and assay-admission settings
envs/                            per-rule conda environments
nanoglobin/
  assay.py                       assay profile validation and in-silico PCR
  cohort.py                      selective-report cohort normalisation
  fastq.py                       strict streaming FASTQ/FASTQ.GZ input
  primer_matching.py             quality-aware terminal primer recognition
  product_evidence.py            compiled-product sequence evidence
  molecules.py                   molecule admission orchestration
  molecule_reporting.py          streaming molecule/product QC summaries
docs/
  IMPLEMENTATION_FOUNDATION.md   cohort-to-caller foundation
  MOLECULE_ADMISSION.md          physical molecule-admission contract
  MSC_VERTICAL_SLICE.md          thesis scope versus post-MSc programme
databases/
  naming_layer.csv               IthaGenes-primary catalogue
  variants.csv                   consumer-facing catalogue
  cnvs.csv                       IthaCNVs-derived CNV catalogue
scripts/
  compile_assay.py               compile declared products over haplotype sequences
  admit_amplicon_reads.py        classify and assign real FASTQ molecules
  normalise_cohort_reports.py    preserve report provenance and missingness semantics
  detect_deletions_cigar.py      spanning-read deletion evidence
  coverage_profile.py            binned relative product/depth evidence
  detect_cnv.py                  unresolved gain/loss candidates from coverage bins
  identify_sv.py                 IthaCNVs catalogue comparator
  vep_annotate.py                REST annotation with catalogue-guided transcript choice
  comprehensive_report.py        merge exact declared sample outputs
  patient_summary.py             per-sample genotype/finding summary
  sample_report.py               ranked per-sample research report
  annotation/annotate_variants.py  annotation and report generation
  simulate/                      deterministic development fixtures
tests/                           compiler, cohort, reporting and molecule-admission tests
giab_HG002.sh                    regional GIAB HG002 comparison
run_pipeline_hprc.sh             HPRC multi-sample comparison driver
simulation/                      development truth VCFs
```
