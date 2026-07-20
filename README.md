# NanoGlobin

Variant calling for haemoglobinopathies from Oxford Nanopore long-read data —
α/β-thalassaemias, structural haemoglobin variants (Hb S, Hb C, Hb E), and
copy-number rearrangements. Targets HBA1/HBA2 (chr16) and HBB (chr11);
developed for thalassaemia and generalised to the globin loci.

```
FASTQ → minimap2 → sort/index ─┬→ Clair3 (SNV/indel, phased) → filter → region filter → VEP → annotate
                               ├→ Sniffles ┐
                               ├→ CuteSV   ┤→ identify_sv (IthaCNVs catalogue)
                               └→ coverage_profile (binned depth, HBB-normalised)
                                            └→ comprehensive report → patient summary → clinical report
```

## Status: mid-development

Validation is simulation-based (Badread `nanopore2023`). See `simulation/README.md`.
The pipeline runs anywhere via `--use-conda` (see Reproducibility). The legacy
BlueBEAR module path (`env_sv.sh` / `env_clair3.sh`, and the `*_path` keys in
`config.yml`) still works but is no longer required.

## Reproducibility

Every rule runs in its own pinned conda environment (`envs/*.yaml`), created
automatically by Snakemake. There is no module loading and no shell-switching —
the whole pipeline runs from a single command.

​```bash
# one-time: snakemake (>=8) with conda >=24.7.1 in its environment

# one-time: obtain the Clair3 ONT model (not bundled with the conda package).
# NanoGlobin pins Clair3 1.0.4, which uses v1 (TensorFlow) models. Download
# r1041_e82_400bps_sup_v430 — or the model matching your basecaller — from ONT
# rerio (https://github.com/nanoporetech/rerio), place it under models/, and set
# config.yml `clair3_models:` to that path (e.g. models/r1041_e82_400bps_sup_v430).

# run everything
snakemake --use-conda --cores 8
​```

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

## Results so far

All from simulated data with ground truth constructed by
`scripts/simulate/make_haplotype.py` (deletions) and `scripts/simulate/make_snv.py`
(point mutations), which apply a known variant to the reference and emit a truth VCF.

**Specificity** (wild-type null control — reads simulated from the reference, so
every call is by construction a false positive):

| caller | false positives |
|---|---|
| Sniffles (coverage-scaled default) | 0 |
| CuteSV `min_support=3` | 38 |
| CuteSV `min_support=25` | 0 |
| Clair3 | 1 (1bp homopolymer, AF 12%) — removed by depth/AF filter |

**Sensitivity and limit of detection** (het `-α3.7`, Sniffles):

| depth | support | DR | DV | fraction | GT | SVLEN |
|---|---|---|---|---|---|---|
| ~930× | 398 | 265 | 194 | 43% | 0/1 | −3812 |
| ~465× | 194 | — | — | 42% | 0/1 | −3812 |
| ~93× | 35 | 58 | 35 | 38% | 0/1 | −3812 |
| ~46× | 19 | 34 | 19 | 36% | 0/1 | −3812 |
| ~28× | 10 | 18 | 10 | 36% | 0/1 | −3812 |

Size exact at every depth. Support falls 40-fold; **variant fraction is invariant
at ~40%**. A fixed `min_support` count optimal at 930× (25, = 2.7% of depth) is
83% of depth at 30× and silently misses carriers. The threshold must be a
fraction of median depth, not a count.

**Copy number** (`scripts/coverage_profile.py`, 200bp bins, HBB-normalised):

| sample | truth | summary ratio | binned profile |
|---|---|---|---|
| WT_control | αα/αα | 1.008 | flat 1.0 |
| HET_a37 | -α3.7/αα | 0.872 | 1.0 → **0.50** → 1.0 (right side) |
| HOM_a37 | -α3.7/-α3.7 | 0.849 | 1.0 → **0.00** → 1.0 |
| HET_a42 | -α4.2/αα | 0.87 | **0.50** → 1.0 (left side) |
| HET_fullalpha | --/αα | 0.503 | flat **0.50** |

A single mean over the 8kb HBA window cannot distinguish het from hom `-α3.7`
(0.872 vs 0.849). Binned depth gives 0.50 vs 0.00, and the step edges recover the
breakpoints (chr16:173,384–177,187). The dip *position* discriminates `-α3.7`
(right-sided) from `-α4.2` (left-sided) — the two commonest α-deletions.

`--/αα` is flat: every bin at 0.50, no internal contrast. Within-sample
normalisation compares bins against each other and therefore cannot detect it —
this is the published failure mode where full α-cluster deletions are reported as
αα/αα. Normalising against HBB (chr11, undeleted) gives 0.503.

**Phasing** — Clair3 `--enable_phasing`; het variants emitted as `0|1`.
SRR37686273 carries six het variants across 1.3 kb of HBB, all in cis, read
directly off single molecules.

**Variant filtering** — a planted Hb S (`HBB:c.20A>T`, called correctly at DP 991,
AF 0.465) was silently dropped by a `QUAL>20` filter (QUAL 19.88). ONT QUAL
miscalibrates for SNVs: 28 false calls in the sample sat at DP=2 with AF=1.0 (they
would pass an AF filter), while the real variant sat at DP 991. Populations
separate cleanly on **depth**, not QUAL. Filter is now `FILTER=PASS` +
`FORMAT/DP≥10` + `FORMAT/AF≥0.15`; QUAL dropped.

## Naming and classification layer

Variants are named and classified against a merged catalogue
(`databases/variants.csv`, 2,839 entries) built from three sources, keyed on HGVS
(`GENE:c.notation`):

- **ClinVar** — pathogenicity for 2,814 globin variants (~705 pathogenic/likely
  pathogenic). `scripts/build_clinvar.py`.
- **HbVar** — common names (Hb S, IVS I-110, …) for 863 variants; 795 filled onto
  ClinVar entries lacking a name. `scripts/merge_hbvar.py`.
- **Curated catalogue** — 68 variants with cohort frequencies (n=1,066).

Structural variants are named separately against **IthaCNVs**
(`databases/cnvs.csv`, 258 CNVs, GRCh38.p13) with reciprocal-overlap matching and
subtype resolution.

Annotation is **catalogue-guided**: where VEP returns several transcripts for a
variant, the transcript whose HGVS matches a catalogue entry is preferred. This
resolves promoter/5′UTR variants (e.g. `c.-138C>A`) that canonical-only annotation
leaves as `upstream_gene_variant`, confirmed on a planted `c.-138C>A` truth sample.

## Co-inheritance flagging

Patients co-inheriting **causative** HBA and HBB variants are flagged: co-inherited
α-thalassaemia suppresses HbA2, the diagnostic marker for β-thal trait, so an HBB
carrier can screen normal on HPLC. In this cohort, HBB heterozygotes below the
3.5% HbA2 cutoff rose with α-globin dose — 25.5% (no HBA variant) → 36.5% (HBA het)
→ 43.3% (HBA hom/comp het). The flag fires only when both loci carry causative
variants; benign background variants do not trigger it.

## Known issues

**Thresholds**
- `cutesv_min_support` is a fixed count (25, = 2.5% of ~1000×, calibrated on the
  wild-type null control). Should be a fraction of median depth: 25 reads is 2.7%
  of depth at 930× but 83% at 30×, and het variants sit at ~40%, so a fixed count
  silently misses carriers below ~60×.

**Callers**
- CuteSV reports `./.` even with `--genotype` and the index present (DV counted,
  DR never computed). Zygosity is taken from Sniffles.
- Clair3 scans all 47 chunks of chr11+chr16 (~225 Mbp) for ~29 kb of amplicon
  reads. Needs `--bed_fn`.

**Annotation**
- Promoter/5′UTR HGVS is transcript-dependent: the legacy `c.-` numbering maps to
  a non-canonical isoform, so catalogue-guided selection is required. Variants not
  in the catalogue and >~340 bp upstream stay `upstream_gene_variant` (correct —
  they are intergenic, not the named promoter variants).
- `classify_mutation_type` matches HGVS with explicit bases (`c.25_26delAA`),
  which modern VEP output never produces (`c.25_26del`). β⁰ variants will
  silently classify as `Unclassified`.
- Multi-nucleotide events (`c.126_129delCTTT`) are annotated per-variant rather
  than per-haplotype; needs bcftools CSQ or VEP haplosaurus (phasing is available).

**Clinical logic**
- `annotate_structural` handles HBA only; HBB structural variants are never
  annotated. Its hardcoded `sv_lookup` dict predates `identify_sv.py`/`cnvs.csv`
  and should be replaced by them.
- β⁰/β⁺ classification lists are hardcoded in Python; should be a data file.

**Pipeline**
- WGS mode is wired (`mode: wgs` in config switches the reference) but has never
  been run. No WGS dataset tested.
- `clinical_annotation` declares one output but the script writes several.
- `comprehensive_report.py` and `patient_summary.py` take no arguments and glob
  the filesystem, so simulated samples appear in clinical reports.

## Data sources

- **ClinVar** — `variant_summary.txt.gz`, merged by `scripts/build_clinvar.py`.
- **HbVar** (https://globin.bx.psu.edu/hbvar) — common names, tab-separated export,
  merged by `scripts/merge_hbvar.py`.
- **IthaCNVs** (https://www.ithanet.eu/db/ithacnv) — CNV breakpoints, GRCh38.p13.
  Parsed by `scripts/parse_ithacnv.py`. 258/311 entries have usable coordinates;
  51 are recorded as "Information unclear" and 2 have transposed digits
  (ithaID 298, 3964 — reported upstream).
- `databases/variant_lookup.csv` — curated catalogue, 68 variants, with cohort
  frequencies (n=1,066).

## Layout

```
Snakefile                        pipeline definition
config.yml                       paths, thresholds, target regions
env_sv.sh / env_clair3.sh        legacy BlueBEAR module sets (superseded by --use-conda)
envs/                            per-rule conda environments (pinned)
databases/
  variant_lookup.csv             curated SNV/indel catalogue + cohort frequencies
  variants.csv                   merged catalogue (curated + ClinVar + HbVar)
  cnvs.csv                       IthaCNVs-derived CNV catalogue
scripts/
  identify_sv.py                 CNV matching against cnvs.csv
  coverage_profile.py            binned depth, normalised to a reference region
  vep_annotate.py                VEP annotation with catalogue-guided transcript choice
  comprehensive_report.py        merge caller outputs
  patient_summary.py             per-sample genotype summary
  annotation/annotate_variants.py  clinical annotation + report generation
  parse_ithacnv.py               IthaCNVs HTML → cnvs.csv
  build_clinvar.py               ClinVar merge → variants.csv
  merge_hbvar.py                 fill common names from HbVar export
  simulate/make_haplotype.py     apply a deletion to a reference, emit truth VCF
  simulate/make_snv.py           apply a point mutation, emit truth VCF
simulation/                      truth VCFs and validation results
```
