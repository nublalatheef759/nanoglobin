# NanoThal

Variant calling for α/β-thalassaemia from Oxford Nanopore data. Targets HBA1/HBA2
(chr16) and HBB (chr11).

```
FASTQ → minimap2 → sort/index ─┬→ Clair3 (SNV/indel, phased) → filter → region filter → VEP → annotate
                               ├→ Sniffles ┐
                               ├→ CuteSV   ┤→ identify_sv (IthaCNVs catalogue)
                               └→ coverage_profile (binned depth, HBB-normalised)
                                            └→ comprehensive report → patient summary → clinical report
```

## Status: mid-development

Validation is simulation-based (Badread `nanopore2023`). See `simulation/README.md`.
Not reproducible outside BlueBEAR — paths in `config.yml` are absolute and local.

## Results so far

All from simulated data with ground truth constructed by
`scripts/simulate/make_haplotype.py`, which splices published breakpoints out of
hg38 and emits a truth VCF.

**Specificity** (wild-type null control — reads simulated from the reference, so
every call is by construction a false positive):

| caller | false positives |
|---|---|
| Sniffles (coverage-scaled default) | 0 |
| CuteSV `min_support=3` | 38 |
| CuteSV `min_support=25` | 0 |
| Clair3 | 1 (1bp homopolymer, AF 12%, QUAL 6.38 — removed by `min_quality: 20`) |

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
| HET_a37 | -α3.7/αα | 0.872 | 1.0 → **0.50** → 1.0 |
| HOM_a37 | -α3.7/-α3.7 | 0.849 | 1.0 → **0.00** → 1.0 |
| HET_fullalpha | --/αα | 0.503 | flat **0.50** |

A single mean over the 8kb HBA window cannot distinguish het from hom `-α3.7`
(0.872 vs 0.849). Binned depth gives 0.50 vs 0.00, and the step edges recover the
breakpoints (chr16:173,384–177,187).

`--/αα` is flat: every bin at 0.50, no internal contrast. Within-sample
normalisation compares bins against each other and therefore cannot detect it —
this is the published failure mode where full α-cluster deletions are reported as
αα/αα. Normalising against HBB (chr11, undeleted) gives 0.503.

**SV naming** — `databases/cnvs.csv`, 258 entries parsed from IthaCNVs
(GRCh38.p13), each carrying its `ithaID`. Matching requires reciprocal overlap
≥50% in both directions plus size concordance ±10%. Simulated `-α3.7` type III is
correctly resolved to type III (100% reciprocal) rather than type I (91%).

## Known issues

**Thresholds**
- `cutesv_min_support: 3` is known-wrong (38 false positives at ~1000×) and kept
  only until the fraction-of-depth work lands.
- Threshold should be a fraction of median depth; currently a fixed count.

**Callers**
- CuteSV reports `./.` even with `--genotype` and the index present (DV counted,
  DR never computed). Zygosity is taken from Sniffles.
- Clair3 scans all 47 chunks of chr11+chr16 (~225 Mbp) for ~29 kb of amplicon
  reads. Needs `--bed_fn`.

**Annotation**
- VEP REST returns `unknown` for every `upstream_gene_variant` — i.e. exactly
  where the β-thal promoter variants sit (`c.-151C>T`, `c.-138C>A`). Indels
  return `api_error`. VEP CLI or GeneBe is the likely fix.
- `classify_mutation_type` matches HGVS with explicit bases (`c.25_26delAA`),
  which modern VEP output never produces (`c.25_26del`). β⁰ variants will
  silently classify as `Unclassified`.
- Multi-nucleotide events (`c.126_129delCTTT`) are annotated per-variant rather
  than per-haplotype; needs bcftools CSQ or VEP haplosaurus.
- Common names available for 70 variants only (the study catalogue). IthaGenes
  has ~3,549 but offers no bulk export.

**Clinical logic**
- Co-inheritance flag fires on any HBA+HBB variant pair, not only causative ones.
- `annotate_structural` handles HBA only; HBB structural variants are never
  annotated.
- β⁰/β⁺ classification lists are hardcoded in Python; should be a data file.

**Pipeline**
- WGS mode untested: `sniffles_sv`, `cutesv_sv` and `clair3_call` hardcode the
  amplicon reference instead of using `REFERENCE`.
- `clinical_annotation` declares one output but the script writes several.
- `comprehensive_report.py` and `patient_summary.py` take no arguments and glob
  the filesystem, so simulated samples appear in clinical reports.
- `summary_table` is orphaned — nothing depends on it.

**Environment**
- Dependencies are split across conda and BlueBEAR modules. Clair3's `PYTHONPATH`
  shadows conda's pysam, so Clair3 and the SV callers cannot share a shell (hence
  `env_sv.sh` / `env_clair3.sh`). `--use-conda` with per-rule environments is the
  fix and is not yet done.

## Data sources

- **IthaCNVs** (https://www.ithanet.eu/db/ithacnv) — CNV breakpoints, GRCh38.p13.
  Parsed by `scripts/parse_ithacnv.py`. 258/311 entries have usable coordinates;
  51 are recorded as "Information unclear" and 2 have transposed digits
  (ithaID 298, 3964 — reported upstream).
- `databases/variant_lookup.csv` — curated catalogue, 70 variants, with cohort
  frequencies (n=1,066).

## Layout

```
Snakefile                        pipeline definition
config.yml                       paths, thresholds, target regions
env_sv.sh / env_clair3.sh        module sets (mutually exclusive — see above)
databases/
  variant_lookup.csv             curated SNV/indel catalogue + cohort frequencies
  cnvs.csv                       IthaCNVs-derived CNV catalogue
scripts/
  identify_sv.py                 CNV matching against cnvs.csv
  coverage_profile.py            binned depth, normalised to a reference region
  vep_annotate.py                VEP REST annotation
  comprehensive_report.py        merge caller outputs
  patient_summary.py             per-sample genotype summary
  annotation/annotate_variants.py  clinical annotation + report generation
  parse_ithacnv.py               IthaCNVs HTML → cnvs.csv
  simulate/make_haplotype.py     apply a deletion to a reference, emit truth VCF
simulation/                      truth VCFs and validation results
```
