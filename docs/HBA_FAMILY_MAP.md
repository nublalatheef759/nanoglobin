# HBA family-coordinate and copy-marker map

## Purpose

HBA1, HBA2 and unequal-crossover hybrids are homologous sequences, but their
observations cannot be projected safely by assuming equal numeric offsets are
equivalent. Insertions, deletions, conversion tracts and structural junctions shift
those offsets.

The family map provides a generated coordinate layer for:

- HBA1/HBA2/hybrid copy markers;
- projection of variants and read observations between copy representations;
- annotation of compiled amplicon products;
- identifying intervals that are copy-informative or ambiguous; and
- visualizing evidence supporting copy identity and hybrid transitions.

It is not a clinical allele catalogue and does not assign a nearest structural name.

## Coordinate keys

One sequence is selected as the anchor. Every aligned observation receives either:

- `R:<position>` — a one-based anchor coordinate; or
- `I:<anchor-bases-before-insertion>:<index>` — an inserted source base at a
  deterministic position between anchor bases.

For each source sequence, the map retains source position, anchor position,
operation, anchor base and source base. Deletions remain gaps and insertions remain
explicit coordinates.

## Marker table

At every family key, the software aggregates the allele carried by each sequence and
classifies the site as invariant, copy-informative, presence/absence or
multi-allelic. Alleles unique to one candidate sequence are reported, but “unique
in this input catalogue” must not be confused with population uniqueness.

## Population extension

The word *population* only describes where additional candidate sequences may come
from. HPRC haplotypes or locally resolved UAE haplotypes can be appended to the
FASTA and the map regenerated. The algorithm and coordinate contract do not change.

This makes the resource incremental and testable. There is no need to manually
complete an all-population map before using copy markers in the current assay.

## Command

```bash
python scripts/build_hba_family_map.py \
  --sequences catalogues/hba_copy_sequences.fa \
  --anchor-id HBA2_anchor \
  --positions-tsv results/family_map/positions.tsv \
  --markers-tsv results/family_map/markers.tsv \
  --map-json results/family_map/map.json \
  --summary-json results/family_map/summary.json
```

## Limits

The current implementation uses deterministic minimum-edit pairwise projections to
one anchor. Very long repetitive rearrangements may admit multiple equally good
alignments; those cases require anchored unique flanks, a local graph or explicit
candidate-haplotype alignment. The map must preserve that ambiguity rather than
invent a precise copy origin.
