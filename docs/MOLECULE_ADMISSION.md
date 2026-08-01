# Primer-aware molecule admission

`compile_assay.py` establishes which amplicons each chromosome haplotype can
physically generate. `admit_amplicon_reads.py` connects those compiled products
to observed ONT FASTQ molecules.

It is an evidence stage, not a diploid genotype caller.

## Inputs

1. A versioned assay profile containing declared primer sequences and product
   length ranges.
2. `compiled_products.json` produced by `compile_assay.py` from the same assay
   profile and haplotype catalogue.
3. FASTQ or FASTQ.GZ reads for one sample.

```bash
python scripts/admit_amplicon_reads.py \
  --sample SAMPLE \
  --fastq fastq/SAMPLE.fastq.gz \
  --profile assays/huang_2023.yml \
  --compiled-json results/assay/compiled_products.json \
  --molecules-tsv results/assay/molecules/SAMPLE.molecules.tsv \
  --product-counts-tsv results/assay/molecules/SAMPLE.product_counts.tsv \
  --summary-json results/assay/molecules/SAMPLE.summary.json
```

The corresponding Snakemake target is:

```bash
snakemake --use-conda --cores 8 \
  results/assay/molecules/SAMPLE.summary.json
```

The rule is available only when both `assay_profile` and
`haplotype_catalogue` are configured.

## Physical admission model

Each read is examined in both orientations. Primer recognition is restricted to
read ends and separates two concepts that were previously easy to conflate:

- `Primer.max_mismatches` in the assay profile controls binding to a genomic
  haplotype during in-silico PCR;
- molecule-admission mismatch thresholds control ONT sequencing error while
  recognising primer sequence in a FASTQ read.

The engine supports:

- complete primers after residual adapter/tail sequence;
- partial primers when a read starts or ends inside the oligonucleotide;
- IUPAC primer codes;
- quality-aware endpoint scores;
- reverse-complement molecule orientation; and
- exact-seed indexing with exhaustive fallback, so the seed is a performance
  filter rather than a sensitivity gate.

## Molecule states

Every read receives one state:

- `complete`: a declared forward/reverse primer pair and compatible product
  length are observed;
- `one_ended`: only one terminal primer is supported;
- `primer_dimer_candidate`: a declared primer pair is present but the span is
  too short for the product;
- `unexpected_length`: a declared pair produces a span outside its product
  contract;
- `chimera_candidate`: the two terminal primers do not form any declared
  product; or
- `off_target`: no declared terminal primer is supported.

None of the unresolved states is converted to a reference genotype.

## Read-to-compiled-product evidence

A complete oriented molecule is trimmed to its observed primer boundaries and
compared with all compatible compiled product sequences for that product ID.
RapidFuzz supplies a native Levenshtein distance. NanoGlobin emits:

- edit distance and normalized edit rate;
- a homogeneous-error log-likelihood used only for ranking;
- the best product-sequence hash;
- all haplotypes that compile to that exact sequence; and
- the margin to the next distinct sequence.

The assignment state is one of:

- `unique_sequence`;
- `equivalent_haplotypes` when multiple chromosome haplotypes yield the same
  amplicon sequence;
- `low_margin` when distinct compiled sequences fit similarly;
- `poor_sequence_fit` when the catalogue does not explain the read;
- `no_compiled_sequence`; or
- `ambiguous_product_definition` when terminal evidence cannot separate two
  declared product definitions.

This deliberately preserves assay-induced ambiguity. The emitted
`sequence_log_likelihood` is not yet a calibrated ONT pair-HMM, and read counts
are not yet converted into a diploid genotype posterior.

## PCR primer incorporation

The assay compiler now emits the sequence of the amplified molecule, not a raw
slice of the genomic template. Tolerated template differences under a primer are
overwritten by the oligonucleotide sequence during PCR. Haplotypes differing
only under an incorporated primer can therefore become observationally
indistinguishable, which is the physically correct result.

## Outputs

`*.molecules.tsv` retains one inspectable row per read, including primer evidence,
orientation, molecule state, product candidates, sequence compatibility and the
reason for ambiguity or failure. For reverse-oriented reads, trim coordinates are
reported in the oriented molecule coordinate system used for sequence comparison.

`*.product_counts.tsv` summarizes endpoint candidates, assigned reads, sequence
ambiguity, poor fits and observed product lengths. Exact length medians are computed
from bounded integer histograms rather than retaining all read observations in
memory. Counts are assay-product abundances, not genomic allele counts.

`*.summary.json` reports assignment, off-target, one-ended, dimer and chimera
fractions. These quantities are directly usable for exploratory assay QC on
unlabelled public amplicon data.

## Claim boundary

This layer can establish that NanoGlobin recognises the physical products in a
real assay and can falsify incompatible haplotype/product hypotheses. It does
not establish sensitivity, specificity, allele dropout, copy-number accuracy or
clinical concordance without assay-matched independent truth.
