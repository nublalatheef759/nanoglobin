# NanoGlobin

![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A58.0-brightgreen)
![Python](https://img.shields.io/badge/python-3.10-blue)
![License](https://img.shields.io/badge/license-MIT-yellow)

Variant calling for haemoglobinopathies from Oxford Nanopore long-read data —
α/β-thalassaemias, structural haemoglobin variants (Hb S, Hb C, Hb E), and
copy-number rearrangements. Targets HBA1/HBA2 (chr16) and HBB (chr11);
developed for thalassaemia (nanothal) and generalised to the globin loci
(nanoglobin).

```
FASTQ → minimap2 → sort/index ─┬→ Clair3 (SNV/indel, phased) → filter → region filter → VEP ─┐
                               ├→ CIGAR deletion scan (spanning reads, class + zygosity) ────┤
                               ├→ Sniffles ┐                                                 │
                               ├→ CuteSV   ┤→ identify_sv (IthaCNVs catalogue) ──────────────┤
                               ├→ coverage_profile (binned depth) → detect_cnv ──────────────┤
                               │                                                             ↓
                               │                              comprehensive report → sample report
                               │                                                    → clinical annotation
                               └→ bcftools csq (haplotype-aware consequence; terminal, not merged)
```

## Status

Core methods are benchmarked against real long-read data using two independent
truth sources (see Benchmarking). The pipeline runs anywhere via `--use-conda`.

> **Running the pipeline after adding/editing a rule:** Snakemake treats a changed
> rule as a reason to rebuild everything upstream. Always pass
> `--rerun-triggers mtime` so only genuinely stale targets rebuild, and put the
> target *after* a `--` separator (the flag otherwise swallows it):
> `snakemake --use-conda --cores 8 --rerun-triggers mtime -- <target>`.

---

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

# one-time: reference data (not bundled -- ~3.5 GB)
# GRCh38 primary assembly, e.g. from UCSC:
#   wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#   gunzip hg38.fa.gz && samtools faidx hg38.fa
# amplicon-mode reference (chr11 + chr16 only):
#   samtools faidx hg38.fa chr11 chr16 > reference/hg38_globin.fa
#   samtools faidx reference/hg38_globin.fa
# annotation for bcftools csq (Ensembl 110 GFF3, chr11 + chr16, UCSC contig names):
#   wget https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/\
#        Homo_sapiens.GRCh38.110.gff3.gz -O reference/ensembl.110.gff3.gz
#   zcat reference/ensembl.110.gff3.gz \
#     | awk 'BEGIN{OFS="\t"} /^#/{print; next} $1=="11"{$1="chr11"; print} $1=="16"{$1="chr16"; print}' \
#     | bgzip > reference/ensembl_globin.gff3.gz
#   tabix -p gff reference/ensembl_globin.gff3.gz

# place reads as fastq/<sample>.fastq and list them in config.yml `samples:`

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

---

## Benchmarking

Benchmarked on real Oxford Nanopore data against two independent truth sources —
HPRC diploid assemblies and the GIAB HG002 benchmark. DRAGEN calls are used as a
**comparator** for deletion carriers, not as independent truth. No single truth
set carries the whole claim.

**1. SNV/indel — assembly and benchmark truth, GA4GH-standard hap.py.**
Pipeline calls benchmarked using hap.py (jmcdani20/hap.py v0.3.12, the
GA4GH-standard benchmarking tool), pooled across six samples: five HPRC genomes
with diploid assembly-derived truth (HG02071, HG02083, HG02514, HG02074,
HG02622) plus HG002 against GIAB v4.2.1 high-confidence truth. HG002 therefore
contributes to the pooled figures below and is not an additional independent
sample:

| region | type | recall | precision |
|---|---|---|---|
| HBA/HBB core | SNV | 96.6% | 99.3% |
| HBA/HBB core | INDEL | 63.6% | 73.7% |
| β-cluster gene bodies | SNV | 99.2% | 99.2% |
| β-cluster gene bodies | INDEL | 100% | 100% |

SNV calling is strong across both the core loci and the extended β-cluster gene
bodies. Indels are the primary limitation (~63% recall), reflecting known ONT
indel error modes, with false negatives clustering at a recurring low-complexity
tandem-repeat region (chr16:171,206–171,221). Across the full extended interval
(including intergenic and locus-control-region sequence) SNV precision falls to
94.1%, driven by false positives in repetitive/regulatory intergenic and LCR
sequence rather than in the globin genes themselves.

**2. SNV/indel — GIAB HG002, region-restricted concordance check.**
Against the NIST v4.2.1 benchmark, restricted to globin targets intersected with
the GIAB high-confidence BED: **23/23 concordant, 0 FN, 0 FP.** This is a
region-restricted check on 23 small variants, not a sensitivity estimate, and
HG002 is the same sample included in the pooled rows above. HG002 is clean in
the globin region (no calls at the tandem-repeat positions).

**3. Deletions — DRAGEN-concordant carriers (CIGAR method), comparator only.**
α-globin deletion carriers with ONT reads, typed by the CIGAR method (below) and
checked against DRAGEN calls. DRAGEN is a **comparator, not independent truth** —
agreement is concordance between two callers, not proof of correctness. No
false-positive deletion call was made in the four DRAGEN-labelled αα/αα controls;
their individual sample identifiers were not retained in the result snapshot:

| sample | comparator label | detected | zygosity |
|---|---|---|---|
| NA21106 | -α4.2/αα | 4257 bp = -α4.2 | het (0.58) |
| HG00642 | -α3.7/αα | 3804 bp = -α3.7 | het (0.71) |
| HG03136 | -α3.7/-α3.7 | 3804 bp = -α3.7 | hom (1.00) |
| HG00735 | ααα3.7/αα | none called | no confident call (2 reads, below threshold) |
| HG03862 | --/αα | none called | no confident call; deferred to coverage |
| normals ×4 (IDs not retained) | αα/αα | none | — |

Synthetic positives (`scripts/simulate/`) were used during development to build
and unit-test the callers against known-answer cases before real data; they are
development scaffolding, not a benchmarking pillar.

---

## α-globin deletion detection (CIGAR method)

`scripts/detect_deletions_cigar.py`. ONT reads long enough to span an entire
α-deletion carry it as a single large `D` operation in their CIGAR string. The
method:

- **detects** the deletion from the CIGAR of spanning reads (no false-positive
  calls across the four αα/αα controls tested);
- **types** it by breakpoint size against `cnvs.csv` size columns, distinguishing
  the `-α3.7` and `-α4.2` deletion classes (the ~450 bp size difference is robust to
  CIGAR-size variance and is corroborated by independent assembly-derived breakpoint
  sizes). The `-α3.7` subtypes (type I/II/III) differ by ~8 bp and are not resolved
  by breakpoint size; class-level reporting is used, as these subtypes are clinically
  equivalent (single-gene α+-deletion);
- **calls zygosity** from the deletion-read fraction *at the breakpoint*
  (het ~0.5–0.7, hom 1.0);
- **makes no deletion call on a triplication** — `ααα3.7` yields no confident
  CIGAR deletion (HG00735), the correct output. The method does **not** identify
  the gain: a triplication and a normal sample return the same verdict, because
  the extra near-identical α-copy maps to the same coordinates and adds no unique
  depth (23.1x vs 24.1x in a normal control).

This is the long-read advantage short reads cannot structurally achieve: precise
class (-α3.7 vs -α4.2) and zygosity within read length. Deletions exceeding read length
(`--`, ~10.5 kb, spanned by too few reads) defer to the coverage method.

**Coverage method.** Binned depth handles large (`--`) deletions the CIGAR method
cannot span. It exists in **two normalisation frames**, and which one applies
depends on the library design — neither covers both.

*HBB-normalised* (`scripts/coverage_profile.py` → `detect_cnv.py`): HBA bin depth
divided by median HBB depth, called against a fixed diploid baseline of 1.0. HBB
lies on chr11, so the yardstick survives even when the entire α locus is deleted.
This is what recovers `--/αα` (ratio 0.504), the genotype causing Hb Bart's
hydrops fetalis, which within-region normalisation reports as normal. Validated on
simulated genotypes: `-α3.7` het and hom, `-α4.2` het, `--/αα`, and a triplication
(1.80× over the duplicated span) all called correctly with a clean wild-type
control. The fixed baseline assumes HBA/HBB ≈ 1.0, which holds for targeted
libraries but not for genome-wide data.

*Panel- and flank-normalised* (`panel_normalise.py`, `detect_alpha_cnv.py`): bin
depth self-normalised against a stable flanking region (chr16:1,000,000–1,100,000)
and compared with a panel of normal `αα/αα` samples. This is required for real WGS,
where HBA/HBB ranges 0.46–0.98 across normal samples because of GC and capture
differences, so a fixed baseline of 1.0 misfires. A deletion is called only where
the deepest window falls below 0.75× the sample's own flanking baseline, which
prevents a uniformly offset profile being read as a loss. Validated on five HPRC
samples with DRAGEN comparator labels: 4/4 deletion carriers detected (relative
0.00–0.72), 2/2 normal controls clean (0.89, 0.94), and the triplication carrier
HG00735 correctly returning no localised event (0.89). This frame requires
coverage outside the globin loci and cannot run on targeted or simulated libraries.

Because HBA1/HBA2 lie in a segmental duplication, neither frame resolves zygosity
reliably for large deletions (het and hom depth overlap), and neither detects
triplications: an extra near-identical α copy adds no unique depth (HG00735,
23.1× vs 24.1× in a normal control). The CIGAR and coverage methods are
complementary — CIGAR for precise typing within read length, coverage for
deletions beyond it.

---

## Naming and classification layer

Variants are named and classified against an **IthaGenes-primary catalogue**
(`databases/naming_layer.csv` → converted to `databases/variants.csv`, 4,027
entries; 3,676 with coordinates), keyed on HGVS (`GENE:c.notation`). This supersedes the earlier
ClinVar-derived catalogue (2,839 entries) — it is a strict superset (+846
variants) that uses expert thalassaemia curation as the classification spine:

- **IthaGenes (ITHANET)** — 2,423 globin variants (2,367 Causative + 56 Neutral),
  the curation spine. Loci: β 1,168; α2 421; α1 283; α-ambiguous 205; δ 193;
  Gγ 88; Aγ 58; plus multi-gene and hybrid entries.
  `databases/build_naming_layer.py`.
- **ClinVar** — 1,604 additional variants not in IthaGenes, plus classification
  cross-reference where both hold a variant.
- Source per variant: 1,411 IthaGenes-only, 1,012 IthaGenes+ClinVar,
  1,604 ClinVar-only (4,027 total; 3,676 with coordinates).

Classification is **Functionality-first**: `sample_report.py:tier()` checks the
IthaGenes Functionality (`Causative`) before ClinVar, so IthaGenes causative
variants tier as pathogenic even where ClinVar is silent or VUS-heavy — the whole
point of the IthaGenes-primary design. 1,394/1,411 IthaGenes-only variants tier
as causative (tier 1). No ClinVar terms are fabricated for IthaGenes-only
variants; a blank ClinVar column is honest and does not demote them.

**GRCh38 coordinates** were derived via the VariantValidator REST API and
validated against 123 known coordinates (45 general + 78 HBB promoter) with zero
mismatches before being trusted on unknowns. Coverage: **99.8% (3,676/3,685)** of coordinate-derivable entries. A further
351 entries are retained for naming only and carry no coordinate: 342 HGVS-keyed
extensions plus nine documented edge cases (protein-only notation,
boundary-spanning deletions, a transcript numbering discrepancy). Unresolved
states are recorded in `coord_status`, never force-filled.
Derivation paths, each validated against known answers:

- direct GRCh38 VCF block (in-transcript variants);
- genomic re-submit suggestion parsed from the validation warning, for
  out-of-transcript substitutions (promoter `c.-`, 3′UTR `c.*`) — VariantValidator
  does the minus-strand arithmetic;
- two-step transcript→genomic query for out-of-transcript indels;
- component-split for compound `[a;b]` alleles (each component resolved, coordinates
  stored semicolon-joined);
- HGVS cleaning (redundant deleted bases, size suffixes) validated 30/30.

The 9 unresolved entries are documented in a `coord_status` column with specific
reasons (protein-notation with no unique nucleotide coordinate, dual-gene-uncertain
HGVS, large boundary-spanning deletions) — flagged for manual curation, never
force-filled with a guessed coordinate. Two coordinate errors found in IthaCNVs
during this work were reported to and corrected by ITHANET, which also added CSV
export to IthaGenes and IthaCNVs in response.

Structural variants are named against **IthaCNVs** (`databases/cnvs.csv`, 311 CNVs,
GRCh38.p13) with reciprocal-overlap matching at the deletion-class level (-α3.7 vs -α4.2; subtypes not resolved by overlap), for gains as
well as deletions.

Annotation is **catalogue-guided**: where VEP returns several transcripts, the one
whose HGVS matches a catalogue entry is preferred, resolving promoter/5′UTR
variants that canonical-only annotation leaves as `upstream_gene_variant`.
Haplotype-aware consequences are additionally produced by **bcftools csq**
(`rule csq_annotate`, Ensembl GFF3), which reports multi-variant haplotypes on a
single allele rather than per-variant.

---

## Co-inheritance flagging

Patients co-inheriting **causative** HBA and HBB variants are flagged. Co-inheritance
alters red-cell indices and complicates interpretation of screening results, so the
combined genotype is surfaced rather than reported as two independent findings. The
flag fires only when both loci carry causative variants; benign background variants
do not trigger it.

---

## Per-sample clinical report

`results/sample_report.csv` gives one row per sample for high-volume review:
a descriptive Result flag (causative / conflicting / VUS-only / none — a
description of what was found, not a diagnosis), the primary finding's common
name, HGVS, and zygosity, variant counts by significance tier, and a full ranked
variant list. SNVs and structural variants (deletions and catalogue-matched
gains) are unified into one view. Nothing is hidden: benign and VUS variants sort
below causative ones but remain visible — ClinVar "Conflicting" variants are
surfaced explicitly.

---

## Variant filtering and phasing

**Filtering** — a planted Hb S (`HBB:c.20A>T`, called correctly at DP 991,
AF 0.465) was silently dropped by a `QUAL>20` filter (QUAL 19.88). ONT QUAL
miscalibrates for SNVs: 28 false calls sat at DP=2 with AF=1.0 while the real
variant sat at DP 991. Populations separate on **depth**, not QUAL. Filter is now
`FILTER=PASS` + `FORMAT/DP≥10` + `FORMAT/AF≥0.15`; QUAL dropped.

**Phasing** — Clair3 `--enable_phasing`; het variants emitted as `0|1`.
Benchmarked against HPRC assembly-derived truth with `whatshap compare` (v2.8,
`run_whatshap_compare.sh`): 256 phased heterozygous variant pairs assessed
across five samples, **zero switch errors** and zero Hamming distance. chr11
only — too few heterozygous variants in the α-globin truth regions to assess.
SRR37686273 carries six het variants across 1.3 kb of HBB, all in cis, read
directly off single molecules.

---

## Known issues

**Callers**
- CuteSV reports `./.` even with `--genotype` and the index present (DV counted,
  DR never computed). Zygosity is taken from Sniffles.
- Clair3 needs `--bed_fn` to avoid scanning all of chr11+chr16 for amplicon reads.

**Annotation**
- Promoter/5′UTR HGVS is transcript-dependent; catalogue-guided selection is
  required. Variants not in the catalogue and >~340 bp upstream stay
  `upstream_gene_variant` (correct — intergenic, not named promoter variants).
- `vep_annotate.py` queries the Ensembl VEP REST endpoint. Requests are retried
  up to four times with exponential backoff, honouring `Retry-After`; before
  this was added a single transient failure marked a variant `api_error` for the
  whole run, which is why a different set of variants failed on each execution.
  A residual failure mode remains: Ensembl sometimes returns an HTML error page
  with a 200 status, which the retry logic cannot detect, so one to three
  variants per run may still record `api_error`. Affected variants are recorded
  rather than silently dropped. Annotation depends on a live external service
  and is not pinned by the conda environments; caching successful annotations,
  batching requests, or using a local VEP install would remove this dependency.

**Clinical logic**
- `annotate_structural` handles HBA only; HBB structural variants are not
  annotated by it. Its hardcoded `sv_lookup` is superseded by
  `identify_sv.py`/`cnvs.csv`.
- CIGAR deletion calls are matched to `cnvs.csv` by size **and** position: the
  catalogue entry must be on the same chromosome and the observed breakpoint
  must fall within its interval (±5 kb). Size tolerance is ±2%. Events matching
  no entry are reported as `uncatalogued` rather than forced to a named class.
- β⁰/β⁺ classification lists are a data file (`beta_classification.csv`).

**Pipeline**
- WGS mode is wired (`mode: wgs` switches the reference) but carriers are run
  manually in WGS coordinates; no full WGS dataset run end-to-end.

---

## Data sources

- **IthaGenes / IthaCNVs** (https://www.ithanet.eu) — variant curation, common
  names, functionality, and CNV breakpoints (GRCh38.p13). `build_naming_layer.py`,
  `refresh_cnvs.py`, `parse_ithacnv.py`.
- **ClinVar** — `variant_summary.txt.gz`, merged by `scripts/build_clinvar.py`.
- **HbVar** (https://globin.bx.psu.edu/hbvar) — common names, merged by
  `scripts/merge_hbvar.py`.
- **VariantValidator** (https://rest.variantvalidator.org) — GRCh38 coordinate
  derivation (MANE transcripts NM_000518.5 HBB, NM_000558.5 HBA1, NM_000517.6 HBA2).

---

## Layout

```
LICENSE                          MIT
Snakefile                        pipeline definition (incl. csq_annotate, cigar_deletions)
config.yml                       paths, thresholds, target regions
envs/                            per-rule conda environments (pinned)
databases/
  naming_layer.csv               IthaGenes-primary catalogue (4027, coord_status column)
  variants.csv                   consumer-facing catalogue (naming_layer, 18-col schema)
  cnvs.csv                       IthaCNVs-derived CNV catalogue (311)
  build_naming_layer.py          build the IthaGenes-primary catalogue
  refresh_cnvs.py                refresh cnvs.csv from IthaCNVs export
scripts/
  detect_deletions_cigar.py      CIGAR-based α-deletion typing + zygosity
  panel_normalise.py             flanking-region coverage normalisation
  detect_alpha_cnv.py            localised-window α-CNV calling
  coverage_profile.py            binned depth, normalised to a reference region
  detect_cnv.py                  copy-number gain/loss from coverage bins
  identify_sv.py                 CNV matching against cnvs.csv
  vep_annotate.py                VEP annotation with catalogue-guided transcript choice
  comprehensive_report.py        merge caller outputs
  patient_summary.py             per-sample genotype summary
  sample_report.py               ranked per-sample clinical report (Functionality-first tiering)
  annotation/annotate_variants.py  clinical annotation + report generation
  naming_layer_to_variants.py    naming_layer → variants.csv (18-col drop-in)
  fill_coordinates.py            VariantValidator coordinate derivation (all paths)
  validate_variantvalidator.py   validate-before-trust gate (general)
  validate_warning_path.py       validate promoter warning-path (78 knowns)
  validate_indel_cleaning.py     validate indel HGVS cleaning (30 knowns)
  finalize_coordinates.py        flag unresolved coordinates with documented reasons
  fill_remaining.py              boundary-spanning deletion recovery
  build_clinvar.py               ClinVar merge
  merge_hbvar.py                 fill common names from HbVar
  simulate/                      development test-case generators (haplotype/SNV/triplication)
giab_HG002.sh                    GIAB HG002 benchmark
run_pipeline_hprc.sh             HPRC multi-sample benchmarking driver
simulation/                      development truth VCFs
```

