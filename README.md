# NanoGlobin

Variant calling for haemoglobinopathies from Oxford Nanopore long-read data —
α/β-thalassaemias, structural haemoglobin variants (Hb S, Hb C, Hb E), and
copy-number rearrangements. Targets HBA1/HBA2 (chr16) and HBB (chr11);
developed for thalassaemia (nanothal) and generalised to the globin loci
(nanoglobin).

```
FASTQ → minimap2 → sort/index ─┬→ Clair3 (SNV/indel, phased) → filter → region filter → VEP → annotate
                               ├→ bcftools csq (haplotype-aware consequence)
                               ├→ CIGAR deletion scan (α-globin deletions from spanning reads)
                               ├→ Sniffles ┐
                               ├→ CuteSV   ┤→ identify_sv (IthaCNVs catalogue)
                               └→ coverage_profile (binned depth, flanking-normalised)
                                            └→ comprehensive report → patient summary → clinical report
```

## Status

Core methods validated against real long-read data with three independent truth
sources (below). The pipeline runs anywhere via `--use-conda` (see
Reproducibility). The legacy BlueBEAR module path still works but is not required.

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

## Validation

Validated on real Oxford Nanopore data against three independent truth sources —
no single truth set carries the whole claim.

**1. SNV/indel — HPRC assembly truth (5 genomes), GA4GH-standard hap.py.**
Pipeline calls benchmarked against HPRC diploid assemblies using hap.py
(jmcdani20/hap.py v0.3.12, the GA4GH-standard benchmarking tool), pooled across
5 genomes (HG02071, HG02083, HG02514, HG02074, HG02622):

| region | type | recall | precision |
|---|---|---|---|
| HBA/HBB core | SNV | 96.6% | 99.3% |
| HBA/HBB core | INDEL | 63.6% | 73.7% |
| β-cluster gene bodies | SNV | 96.9% | 99.4% |
| β-cluster gene bodies | INDEL | 63.3% | 91.2% |

SNV calling is strong across both the core loci and the extended β-cluster gene
bodies. Indels are the primary limitation (~63% recall), reflecting known ONT
indel error modes, with false negatives clustering at a recurring low-complexity
tandem-repeat region (chr16:171,206–171,221). Across the full extended interval
(including intergenic and locus-control-region sequence) SNV precision falls to
94.1%, driven by false positives in repetitive/regulatory intergenic and LCR
sequence rather than in the globin genes themselves.

**2. SNV/indel — GIAB gold standard (HG002).**
Against NIST v4.2.1 benchmark, restricted to globin high-confidence regions:
**23/23 true positives, 0 FN, 0 FP = 100% / 100%.** HG002 is clean in the globin
region (no calls at the tandem-repeat positions), complementing the HPRC breadth.

**3. Deletions — DRAGEN-concordant carriers (CIGAR method).**
α-globin deletion carriers with ONT reads, typed by the CIGAR method (below) and
checked against DRAGEN calls — 0 false positives:

| sample | truth | detected | zygosity |
|---|---|---|---|
| NA21106 | -α4.2/αα | 4257 bp = -α4.2 | het (0.58) |
| HG00642 | -α3.7/αα | 3804 bp = -α3.7 | het |
| HG03136 | -α3.7/-α3.7 | 3804 bp = -α3.7 | hom (1.00) |
| HG00735 | ααα3.7/αα | (triplication) | correctly excluded |
| HG03862 | --/αα | 10553 bp | deferred to coverage (below read span) |
| normals ×4 | αα/αα | none | — |

Synthetic positives (`scripts/simulate/`) were used during development to build
and unit-test the callers against known-answer cases before real data; they are
development scaffolding, not a validation pillar.

## α-globin deletion detection (CIGAR method)

`scripts/detect_deletions_cigar.py`. ONT reads long enough to span an entire
α-deletion carry it as a single large `D` operation in their CIGAR string. The
method:

- **detects** the deletion from the CIGAR of spanning reads (0 false positives);
- **types** it by breakpoint size against `cnvs.csv` size columns, distinguishing
  the `-α3.7` and `-α4.2` deletion classes (the ~450 bp size difference is robust to
  CIGAR-size variance and is corroborated by independent assembly-derived breakpoint
  sizes). The `-α3.7` subtypes (type I/II/III) differ by ~8 bp and are not resolved
  by breakpoint size; class-level reporting is used, as these subtypes are clinically
  equivalent (single-gene α+-deletion);
- **calls zygosity** from the deletion-read fraction *at the breakpoint*
  (het ~0.5–0.7, hom 1.0);
- **excludes triplications** — an extra near-identical α-copy maps to the same
  coordinates and produces no depth change, so `ααα3.7` is correctly *not* called
  as a deletion (confirmed on HG00735).

This is the long-read advantage short reads cannot structurally achieve: precise
class (-α3.7 vs -α4.2) and zygosity within read length. Deletions exceeding read length
(`--`, ~10.5 kb, spanned by too few reads) defer to the coverage method.

**Coverage method** (`scripts/coverage_profile.py`, `panel_normalise.py`,
`detect_alpha_cnv.py`): binned depth normalised against a *stable flanking region*
(chr16:1,000,000–1,100,000), measuring the deepest localised window. This handles
large (`--`) deletions the CIGAR method cannot span, but — because HBA1/HBA2 lie
in a segmental duplication — coverage cannot reliably resolve zygosity (het/hom
depth overlap) and misclassifies triplications as deletions. The two methods are
complementary: CIGAR for precise typing within read length, flanking-normalised
coverage for large deletions.

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
variants tier as pathogenic even where ClinVar is silent or VUS-heavy — the whole
point of the IthaGenes-primary design. 1,055/1,069 IthaGenes-only variants tier
as causative (tier 1). No ClinVar terms are fabricated for IthaGenes-only
variants; a blank ClinVar column is honest and does not demote them.

**GRCh38 coordinates** were derived via the VariantValidator REST API and
validated against 123 known coordinates (45 general + 78 HBB promoter) with zero
mismatches before being trusted on unknowns. Coverage: **99.8% (3,676/3,685)**.
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

## Co-inheritance flagging

Patients co-inheriting **causative** HBA and HBB variants are flagged: co-inherited
α-thalassaemia suppresses HbA2, the diagnostic marker for β-thal trait, so an HBB
carrier can screen normal on HPLC. In this cohort, HBB heterozygotes below the
3.5% HbA2 cutoff rose with α-globin dose — 25.5% (no HBA variant) → 36.5% (HBA het)
→ 43.3% (HBA hom/comp het). The flag fires only when both loci carry causative
variants; benign background variants do not trigger it.

## Per-sample clinical report

`results/sample_report.csv` gives one row per sample for high-volume review:
a descriptive Result flag (causative / conflicting / VUS-only / none — a
description of what was found, not a diagnosis), the primary finding's common
name, HGVS, and zygosity, variant counts by significance tier, and a full ranked
variant list. SNVs and structural variants (deletions and catalogue-matched
gains) are unified into one view. Nothing is hidden: benign and VUS variants sort
below causative ones but remain visible — ClinVar "Conflicting" variants are
surfaced explicitly.

## Variant filtering and phasing

**Filtering** — a planted Hb S (`HBB:c.20A>T`, called correctly at DP 991,
AF 0.465) was silently dropped by a `QUAL>20` filter (QUAL 19.88). ONT QUAL
miscalibrates for SNVs: 28 false calls sat at DP=2 with AF=1.0 while the real
variant sat at DP 991. Populations separate on **depth**, not QUAL. Filter is now
`FILTER=PASS` + `FORMAT/DP≥10` + `FORMAT/AF≥0.15`; QUAL dropped.

**Phasing** — Clair3 `--enable_phasing`; het variants emitted as `0|1`.
SRR37686273 carries six het variants across 1.3 kb of HBB, all in cis, read
directly off single molecules.

## Known issues

**Callers**
- CuteSV reports `./.` even with `--genotype` and the index present (DV counted,
  DR never computed). Zygosity is taken from Sniffles.
- Clair3 needs `--bed_fn` to avoid scanning all of chr11+chr16 for amplicon reads.

**Annotation**
- Promoter/5′UTR HGVS is transcript-dependent; catalogue-guided selection is
  required. Variants not in the catalogue and >~340 bp upstream stay
  `upstream_gene_variant` (correct — intergenic, not named promoter variants).
- `vep_annotate.py` annotates insertions via the Ensembl VEP region endpoint,
  which can return an API error for some indel formats. Not triggered by any
  variant in the current test set.

**Clinical logic**
- `annotate_structural` handles HBA only; HBB structural variants are not
  annotated by it. Its hardcoded `sv_lookup` is superseded by
  `identify_sv.py`/`cnvs.csv`.
- β⁰/β⁺ classification lists are a data file (`beta_classification.csv`).

**Pipeline**
- WGS mode is wired (`mode: wgs` switches the reference) but carriers are run
  manually in WGS coordinates; no full WGS dataset run end-to-end.
- `comprehensive_report.py` and `patient_summary.py` take no arguments and glob
  the filesystem, so simulated samples can appear in clinical reports.

## Data sources

- **IthaGenes / IthaCNVs** (https://www.ithanet.eu) — variant curation, common
  names, functionality, and CNV breakpoints (GRCh38.p13). `build_naming_layer.py`,
  `refresh_cnvs.py`, `parse_ithacnv.py`.
- **ClinVar** — `variant_summary.txt.gz`, merged by `scripts/build_clinvar.py`.
- **HbVar** (https://globin.bx.psu.edu/hbvar) — common names, merged by
  `scripts/merge_hbvar.py`.
- **VariantValidator** (https://rest.variantvalidator.org) — GRCh38 coordinate
  derivation (MANE transcripts NM_000518.5 HBB, NM_000558.5 HBA1, NM_000517.6 HBA2).

## Layout

```
Snakefile                        pipeline definition (incl. csq_annotate, cigar_deletions)
config.yml                       paths, thresholds, target regions
envs/                            per-rule conda environments (pinned)
databases/
  naming_layer.csv               IthaGenes-primary catalogue (3685, coord_status column)
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
giab_HG002.sh                    GIAB HG002 gold-standard validation
run_pipeline_hprc.sh             HPRC multi-sample validation driver
simulation/                      development truth VCFs
```

