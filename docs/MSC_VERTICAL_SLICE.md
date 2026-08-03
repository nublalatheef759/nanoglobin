# Current research programme, thesis outputs and validation gates

The previous version of this document incorrectly pushed the HBA family-coordinate
map and chromosome-haplotype posterior into a generic “post-MSc programme.” That
framing was too conservative and technically misleading. Both are ordinary,
important parts of the current method and can be implemented incrementally from
well-established alignment, paralog-phasing, mixture-likelihood and amplicon-assay
precedents.

The useful distinction is not **MSc versus post-MSc**. It is:

1. what is implemented and executable;
2. what is active model engineering and can produce thesis results; and
3. what requires new biological data or governance before performance or deployment
   claims are justified.

## The current computational method

```text
candidate chromosome-haplotype sequences
    -> declared primer/product assay profile
    -> in-silico PCR and assay-observability matrix
    -> HBA family-coordinate and copy-marker projection
    -> primer-aware ONT molecule admission
    -> product/sequence/artifact evidence
    -> posterior over pairs of chromosome haplotypes
    -> resolved / assay-equivalent / posterior-ambiguous / no-call
```

### Variable copy number is already inside the candidate haplotypes

Humans remain diploid at the chromosome level. An HBA chromosome haplotype may
contain one, two, three or another number of alpha-like copies, with a declared copy
order, hybrid sequence and linked variants. A posterior over pairs of those
chromosome haplotypes is therefore already a variable-copy HBA posterior. It does
not require a separate variable-ploidy model.

The candidate catalogue can begin with normal, common deletion, gain and hybrid
haplotypes and expand as additional HPRC, cohort or orthogonally resolved sequences
become available.

## What the HBA family-coordinate map means

“Full HBA1/HBA2 population equivalence map” was a poor phrase. The implemented
object is an **HBA family-coordinate and marker map**:

- choose one sequence as a coordinate anchor;
- align HBA1, HBA2, hybrid and candidate-copy sequences to it;
- assign homologous bases to shared reference-position keys;
- assign source insertions to explicit insertion keys;
- retain deletions and substitutions rather than forcing equal offsets;
- identify copy-informative, presence/absence and multi-allelic markers; and
- project source-sequence positions and assay products through the same coordinate
  system.

The map is generated from FASTA sequences. Adding population haplotypes means adding
sequences and regenerating the outputs; it is not a manually curated, all-human
prerequisite. Exact compiled product sequences remain the primary read-likelihood
objects. The family map supports interpretation, copy-marker evidence, variant
projection and visualization.

See [`HBA_FAMILY_MAP.md`](HBA_FAMILY_MAP.md).

## Active thesis and method work

The following are current work, not deferred research aspirations:

- expanding the sequence-resolved chromosome-haplotype catalogue;
- generating the HBA family-coordinate and copy-marker map;
- scoring pairs of chromosome haplotypes from molecule evidence;
- modelling product multiplicity and configurable product efficiency;
- explicit dropout and artifact components;
- collapsing assay-equivalent genotype labels before assigning priors;
- returning resolved, ambiguity and no-call states;
- testing synthetic known genotypes and deliberate indistinguishability;
- describing real public amplicon endpoint, length, product and artifact
  distributions; and
- linking the resulting figures and tables directly into the thesis workspace.

The initial likelihood can be transparent and auditable. A pair-HMM, richer ONT
context model or learned likelihood may improve it, but is not a prerequisite for
having a real probabilistic genotyper.

## Calibration is a continuum, not a reason to omit the model

Product efficiency and dropout can advance in layers:

1. declared neutral defaults and sensitivity analysis;
2. estimates from public unlabelled runs for length, artifact and batch dispersion;
3. truth-matched controls for genotype-conditioned efficiencies and dropout;
4. hierarchical lot/laboratory calibration once enough replicated data exist.

The current model exposes every parameter and labels the posterior as uncalibrated
until the appropriate truth-matched controls are available. This is stronger than
hard-coding a support threshold and calling it a genotype model.

## Data- and governance-dependent gates

The following are not postponed because they are computationally difficult. They
require data or permissions that the repository does not currently possess:

- exact AmplideX primer/product or legally accessible assay-matched raw data;
- certified or event-appropriate orthogonal truth for analytical sensitivity and
  specificity;
- blinded held-out clinical evaluation;
- lot, laboratory and operator replication;
- clinical report validation, regulatory interpretation and deployment claims; and
- public release of any resource derived from governed patient data.

Automatic consensus or assembly of coherent unexplained products can also be added
now, but a novel allele name or clinical assertion still requires independent
confirmation.

## Thesis outputs

The thesis should lead with executable results rather than a feature checklist:

- cohort and report data audit;
- canonical reported genotype/haplotype landscape;
- HBA × HBB phenotype analysis where eligible;
- database and nomenclature representation;
- compiled assay coverage and indistinguishability;
- primer-aware molecule QC on real amplicon reads;
- candidate genotype posteriors and explicit failure states;
- synthetic falsification and truth-scoped benchmarking; and
- a transparent account of which claims await assay-matched ground truth.

The Quarto/Typst workspace under `thesis/` renders the same draft to reviewable HTML
and submission-oriented PDF and consumes pipeline outputs directly.
