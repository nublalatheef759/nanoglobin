# Simulation-based validation

All samples generated with Badread v0.4.2, --error_model nanopore2023.
Haplotypes built by scripts/simulate/make_haplotype.py from reference/sim_target.fa
(chr16:165000-182000 + chr11:5220000-5232000, cut from hg38_globin.fa).

## Samples

WT_control      wild-type null control, seed 42, ~1000x
HOM_a37         hom -a3.7,  del chr16:173708-177519 (3812bp), seed 42
HET_a37         het -a3.7,  seeds 1/2, 500x each
HET_{500,100,50,30}x   het -a3.7 coverage series
HET_fullalpha   het --/aa, del chr16:170000-179000 (9001bp), seeds 1/2

## Results

### Copy-number normalisation (HBA vs HBB)
sample          HBA      HBB      ratio   truth
WT_control      969.369  963.662  1.006   aa/aa
HET_fullalpha   490.041  972.640  0.504   --/aa

Within-sample normalisation (Asuragen's published method) reads a flat ~490 across
the whole HBA window on HET_fullalpha and reports aa/aa -- a false negative on the
genotype causing Hb Bart's hydrops fetalis. HBB-normalised depth gives 0.50.

### SV detection (Sniffles)
depth   support  DR   DV   fraction  GT   SVLEN
930x    398      265* 194* ~43%      0/1  -3812
465x    194      -    -    ~42%      0/1  -3812
93x     35       58   35   ~38%      0/1  -3812
46x     19       34   19   ~36%      0/1  -3812
28x     10       18   10   ~36%      0/1  -3812
Size exact at every depth. Fraction invariant; count falls 40-fold.

### Specificity (WT_control)
CuteSV min_support:  3 -> 38 FP | 10 -> 1 | 25 -> 0 | 50 -> 0 | 100 -> 0
Sniffles: 0 FP throughout.
Clair3: 1 call (1bp homopolymer del, AF 12%, QUAL 6.38), removed by min_quality 20.

## Simulating copy-number gains (triplications)

A deletion is simulated by removing a segment (make_haplotype.py). A gain is the
inverse, but must NOT be simulated by running Badread on a tandem-duplicated
reference at fixed Nx -- Badread spreads coverage across the longer template, so
the duplicated span only reaches ~1.3x instead of 2x.

Correct approach: baseline reads from the normal reference, PLUS an extra copy's
worth of reads from just the duplicated segment, so the segment gets genuine 2x:

    # baseline (normal 2-copy depth)
    badread simulate --reference reference/sim_target.fa --quantity 1000x \
        --seed 10 --error_model nanopore2023 | sed 's/^@/@base_/' > base.fq
    # extra copy over the duplicated segment only
    samtools faidx reference/sim_target.fa chr16:173384-177187 > dupseg.fa
    badread simulate --reference dupseg.fa --quantity 1000x \
        --seed 11 --error_model nanopore2023 | sed 's/^@/@extra_/' > extra.fq
    cat base.fq extra.fq > fastq/HOM_triple.fastq

Validated: the binned profile shows ~2x (1.80 vs 0.91 flanking) over the
duplicated span, with sharp boundaries. make_triplication.py emits the DUP truth
VCF for benchmarking.

## Simulating copy-number gains (triplications)
Do NOT run Badread on a tandem-duplicated reference at fixed Nx (coverage spreads
across the longer template; the dup span only reaches ~1.3x). Instead: baseline
reads from the normal reference PLUS an extra copy's worth from just the
duplicated segment, giving a genuine 2x. See make_triplication.py for the truth VCF.
