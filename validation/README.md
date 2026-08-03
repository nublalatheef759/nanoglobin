# Validation data programme

Ground truth is not one VCF and it is not another caller.  NanoGlobin records
truth per sample, chromosome haplotype and event class, together with the assay
that established it.  The resource registry is in [`resources.tsv`](resources.tsv).

## What can be done immediately

### 1. Public assay mechanics

Use `PRJNA1439314` for questions that do not require genotype truth:

- recurrent 5' and 3' endpoints;
- inferred primer/product pairs;
- complete, one-ended, off-target and chimera fractions;
- observed product-length distributions;
- product imbalance and batch dispersion;
- read-assignment residuals; and
- runtime and memory benchmarks.

It must not be used to estimate genotype sensitivity, specificity or calibrated
dropout because no independently established sample-level genotypes or complete
assay contract have been identified.

### 2. Public sequence truth and candidate discovery

Use HPRC raw reads and diploid assemblies to expand and audit the candidate
chromosome-haplotype catalogue.  Each HBA locus must first pass a locus-resolution
audit: unique flanks on both sides, expected copy count and order, junction
sequence, and independent molecule support.  A sequence-continuous contig across a
deletion is expected and is not evidence that the deletion was missed.

Use GIAB HG002 v5.0q only after intersecting the exact HBA/HBB regions and variant
classes with the released benchmark BED.  Report only bases and event classes that
are inside the benchmark.  These WGS resources validate sequence/caller behaviour,
not multiplex-PCR efficiency or allele dropout.

### 3. Executable published assay

[`../assays/huang_2023_ont_long_pcr.yml`](../assays/huang_2023_ont_long_pcr.yml)
transcribes the published primer table from Huang et al.  Compile it against normal,
deletion, gain and hybrid chromosome haplotypes.  This gives a real public assay
contract for:

- expected-product tests;
- assay coverage matrices;
- genotype indistinguishability;
- exact-product simulation; and
- design sensitivity analyses.

The published profile is not AmplideX and does not replace raw reads from the 158
pregenotyped samples.

## Data requests to send now

### Huang et al. 2023

Request:

1. demultiplexed FASTQ for all 158 blinded samples;
2. sample-to-genotype truth table, retaining compound genotypes and normal samples;
3. run, flow-cell, barcode, pool, batch and replicate identifiers;
4. per-product primer/pool assignment and any excluded/no-call reads;
5. truth method for each event: GAP-PCR, RDB, Sanger or MLPA; and
6. permission to release derived aggregate calibration parameters and benchmark
   results.

This is the highest-value near-term request because the paper publishes the exact
primer design and used 158 pregenotyped samples, including deletion, multicopy,
hybrid and HBB classes.

### Yang et al. 2026

Request raw reads, indexed-primer sequences, sample truth, and the validation versus
application cohort split for the 455-sample one-step library-preparation study.
This is a large external challenge set, but its Qitan platform and assay should be
modelled as a separate assay profile.

### Wei et al. 2025

Submit an early-access request for `PRJCA038473` / `HRA011058`, specifically the 32
national reference materials and four MLPA-characterised triplication samples.
The repository has been reported as unavailable until 2027-04-09; verify the live
record and access terms rather than assuming that deposited means downloadable.
The project must not depend on early access being granted.

## Orderable physical controls

A first exact-assay control run can be built from available Coriell/CDC material:

- `CD00027`: heterozygous `--SEA/alpha alpha`, molecularly characterised before
  submission to the CDC repository;
- `NA10799`: independent `--SEA/alpha alpha` family material;
- `NA10798`: heterozygous `--FIL/alpha alpha`;
- `NA10797`: compound `--FIL/--SEA`;
- `NA07407`: compound HBB splice variants, IVS-I-6 T>C and IVS-I-1 G>A; and
- several CDC panel negatives for the specifically tested SEA/FIL and HbS/HbC
  mutations.

Before freezing truth, confirm the current product form, lot, quantity, certificate
and consent/use constraints, then reverify the relevant event with an orthogonal
method.  A panel-negative sample is not automatically a completely normal globin
haplotype.

The largest current physical-control gaps are common `-alpha3.7`, `-alpha4.2`,
anti-3.7/anti-4.2 gains, HKalphaalpha, non-deletional HBA alleles, and beta-cluster
deletions.  Fill these through the Huang/Wei/Yang requests, the laboratory's
historical bank, or newly characterised material.

## Exact AmplideX evaluation

The decisive experiment is an exact-assay, blinded set containing:

- raw AmplideX FASTQ;
- the assay and software version;
- the historical orthogonal result that existed before ONT testing;
- sample-level indication of repeated or related samples;
- all failed, negative, ambiguous and no-call runs;
- run, barcode, pool, operator, lot and date metadata; and
- vendor Reporter output retained as a comparator only.

Preferred truth by event class:

| Event class | Minimum independent evidence |
|---|---|
| HBA/HBB SNV or short indel | bidirectional Sanger or independently validated reference genotype |
| `-alpha3.7`, `-alpha4.2`, `--SEA`, `--FIL` | event-specific GAP-PCR plus junction sequencing |
| copy gain or triplication | MLPA or ddPCR plus junction/copy-order evidence |
| large beta-cluster deletion | MLPA/dosage plus breakpoint-spanning sequence |
| phase or compound genotype | one molecule spanning informative sites, segregation, or orthogonal long-range assay |
| novel product | independent PCR and Sanger/long-read confirmation of the assembled junction |

## Minimum useful experimental design

Do not spend all controls on more copies of the same heterozygote.  The first panel
should deliberately include:

1. multiple independent normal/negative materials;
2. heterozygous and homozygous/compound deletion states;
3. the common alpha-plus and alpha-zero deletion classes;
4. gains/triplications and hybrid alleles;
5. HBA non-deletional SNVs;
6. HBB promoter, splice, coding SNV and indel classes;
7. at least one beta-cluster structural event;
8. deliberate assay-equivalent genotype pairs;
9. serial DNA input and allele-fraction mixtures for limit-of-detection work; and
10. within-run, between-run, lot and operator replicates.

Calibration and evaluation must be separated.  Fit product efficiencies and dropout
on a training subset, freeze the parameters and candidate catalogue, then evaluate
held-out samples blindly.  Report genotype-class confusion, sensitivity,
specificity, no-call rate, ambiguity rate, calibration, coverage by event class,
and every discordance after adjudication.

## Control manifest

The calibration CLI consumes a tab-separated manifest:

```text
sample_id  molecules_tsv  haplotype_1  haplotype_2  stratum  truth_status  truth_source  assay_scope
```

Allowed `truth_status` values are:

- `certified_reference`;
- `orthogonally_resolved`;
- `assembly_resolved`;
- `synthetic`; and
- `comparator_only`.

Allowed assay scopes include `assay_matched`, `synthetic_assay`, `wgs`,
`caller_comparator` and `not_assay_matched`.  Comparator-only and non-assay-matched
controls are excluded from calibration by default.
