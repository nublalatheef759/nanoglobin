#!/usr/bin/env python3
"""
detect_deletions_cigar.py -- long-read deletion detection for the globin regions.

Oxford Nanopore reads are long enough to span whole alpha-globin deletions, so a
deletion appears DIRECTLY inside a single read's alignment as a large CIGAR 'D'
operation, rather than only as a coverage dip. This is the ONT long-read
advantage: short-read methods must infer deletions from depth or split-reads,
whereas a long read shows the deletion breakpoint-to-breakpoint in one alignment.

This detector, validated against DRAGEN-genotyped carriers:
  * scans reads in the target regions for CIGAR deletions above a size threshold
  * clusters them by size and reports the consensus deletion size
  * maps the size to a known deletion type via databases/cnvs.csv (size_min/max)
  * estimates zygosity from the fraction of spanning reads carrying the deletion
    (~0.5 -> heterozygous, ~1.0 -> homozygous)
  * reports NO deletion for normals and for triplications (a gain has no deletion
    signature -- this correctly distinguishes aaa3.7 triplications, which coverage
    methods misclassify)

VALIDATION (5 DRAGEN-genotyped ONT carriers + normal controls):
  -a4.2/aa      -> 4257 bp deletion, fraction 0.57  (het)      correct
  -a3.7/aa      -> 3804 bp deletion, fraction 0.71  (het)      correct
  -a3.7/-a3.7   -> 3804 bp deletion, fraction 1.00  (hom)      correct
  --/aa         -> 10553 bp deletion (see LIMITATION)          correct type
  aaa3.7/aa     -> no deletion (triplication)                  correct
  normals       -> no deletion                                 correct

LIMITATION: a read must be long enough to span the ENTIRE deletion to show it as
one CIGAR 'D'. Deletions up to ~read length (the ~3.8-4.2 kb single-gene
deletions) are captured with good depth; very large deletions (e.g. the ~10.5 kb
-- double-gene deletion) are spanned by fewer reads, so the deletion-read
fraction under-counts and zygosity from this method alone is unreliable for them.
For large deletions, combine with the (flanking-normalised, localised) coverage
method, which detects them robustly by depth. The two methods are complementary:
CIGAR-deletion gives precise size/type + zygosity for deletions within read
length; coverage catches the large ones.
"""
import argparse
import csv
import os
import re
import statistics as stats
import subprocess
import sys


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar_dels(cigar):
    """Yield (ref_offset_at_del_start, del_len) for each D op, tracking ref pos."""
    ref = 0
    for num, op in CIGAR_RE.findall(cigar):
        n = int(num)
        if op in "MDN=X":
            if op == "D":
                yield ref, n
            ref += n
        # I,S,H,P consume no reference


def scan_region(samtools, bam, region, min_del, mapq):
    """Return (dels, reads) where dels = list of (breakpoint_pos, del_len) and
    reads = list of (start, end) for all aligned reads in the region (used later
    to count how many span a given breakpoint)."""
    proc = subprocess.run(
        [samtools, "view", "-q", str(mapq), bam, region],
        capture_output=True, text=True)
    dels = []
    reads = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        pos = int(f[3])
        cigar = f[5]
        if cigar == "*":
            continue
        reflen = sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X")
        read_end = pos + reflen
        reads.append((pos, read_end))
        for off, dlen in parse_cigar_dels(cigar):
            if dlen >= min_del:
                dels.append((pos + off, dlen))
    return dels, reads


def count_spanning(reads, breakpoint_pos, flank=100):
    """How many reads truly cross the breakpoint (start before, end after)."""
    lo = breakpoint_pos - flank
    hi = breakpoint_pos + flank
    return sum(1 for (s, e) in reads if s < lo and e > hi)


def load_cnv_types(cnvs_csv):
    """Return list of (name, size_min, size_max, genes) for typing by size."""
    types = []
    if not os.path.exists(cnvs_csv):
        return types
    with open(cnvs_csv) as fh:
        for r in csv.DictReader(fh):
            try:
                smin = int(float(r.get("size_min", "") or 0))
                smax = int(float(r.get("size_max", "") or 0))
            except ValueError:
                continue
            if smin > 0:
                types.append((r.get("name", "?"), smin, smax or smin,
                              r.get("genes", "")))
    return types


def type_by_size(size, types, tol=0.10):
    """Return CNV records whose size matches, sorted by closeness (best first)."""
    hits = []
    for name, smin, smax, genes in types:
        lo = min(smin, smax) * (1 - tol)
        hi = max(smin, smax) * (1 + tol)
        if lo <= size <= hi:
            # distance from the nearest catalogued size bound
            centre = (smin + smax) / 2
            dist = abs(size - centre)
            hits.append((dist, name, smin, smax, genes))
    hits.sort(key=lambda x: x[0])
    return [(n, a, b, g) for (_, n, a, b, g) in hits]


def write_tsv(path, name, detected, size, mtype, zyg, frac, genotype):
    """Write a single clean result row (with header) for pipeline consumption."""
    if not path:
        return
    header = ("sample\tdeletion_detected\tsize_bp\tmatched_type\t"
              "zygosity\tfraction\tdragen_genotype\n")
    row = (f"{name}\t{detected}\t{size}\t{mtype}\t{zyg}\t{frac}\t{genotype}\n")
    with open(path, "w") as fh:
        fh.write(header)
        fh.write(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--name", default="sample")
    ap.add_argument("--genotype", default="")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--cnvs", default="databases/cnvs.csv")
    ap.add_argument("--regions", default="chr16:170000-180000,chr11:5220000-5240000",
                    help="comma-separated regions to scan")
    ap.add_argument("--min-del", type=int, default=2000,
                    help="minimum CIGAR deletion size to report (bp)")
    ap.add_argument("--min-reads", type=int, default=3,
                    help="minimum supporting reads to call a deletion "
                         "(filters 1-2 read noise; large deletions spanned by "
                         "few reads fall below this and defer to coverage)")
    ap.add_argument("--mapq", type=int, default=20)
    ap.add_argument("--out", default="",
                    help="optional path to write a clean TSV result row "
                         "(sample, detected, size_bp, matched_type, zygosity, "
                         "fraction, genotype)")
    args = ap.parse_args()

    types = load_cnv_types(args.cnvs)
    print(f"# {args.name}" + (f"  (DRAGEN: {args.genotype})" if args.genotype else ""))

    all_dels = []
    all_reads = []
    for region in args.regions.split(","):
        dels, reads = scan_region(args.samtools, args.bam, region,
                                  args.min_del, args.mapq)
        all_dels.extend(dels)
        all_reads.extend(reads)

    # require a minimum number of supporting reads (noise filter)
    if len(all_dels) < args.min_reads:
        if all_dels:
            sizes = sorted(d[1] for d in all_dels)
            print(f"#   {len(all_dels)} read(s) with a large deletion "
                  f"(~{int(stats.median(sizes))}bp) -- below the {args.min_reads}-read "
                  f"threshold.")
            print(f"#   RESULT: no confident CIGAR deletion. Either a normal/gain, "
                  f"or a large deletion spanned by too few reads (confirm by coverage).")
        else:
            print("#   RESULT: no large deletion detected "
                  "(normal, or a gain/triplication -- not depth-distinguishable here)")
        print(f"#   VERDICT: NO CONFIDENT DELETION ({args.name})", file=sys.stderr)
        write_tsv(args.out, args.name, "no", "NA", "none", "NA", "NA", args.genotype)
        return

    sizes = sorted(d[1] for d in all_dels)
    consensus = int(stats.median(sizes))
    n_del = len(all_dels)

    # spanning reads counted at the actual breakpoint (median deletion start)
    breakpoint_pos = int(stats.median(d[0] for d in all_dels))
    spanning = count_spanning(all_reads, breakpoint_pos)
    frac = n_del / spanning if spanning else float("nan")

    if frac >= 0.85:
        zyg = "HOMOZYGOUS (fraction ~1.0)"
    elif frac >= 0.30:
        zyg = "HETEROZYGOUS (fraction ~0.5)"
    else:
        zyg = ("UNCERTAIN (low fraction -- large deletion exceeding typical read "
               "length; confirm zygosity by coverage)")

    hits = type_by_size(consensus, types)

    print(f"#   consensus deletion size: {consensus} bp "
          f"(range {sizes[0]}-{sizes[-1]}, n={n_del} reads)")
    print(f"#   breakpoint ~chr16:{breakpoint_pos}   spanning reads: {spanning}   "
          f"deletion-read fraction: {frac:.2f}")
    if hits:
        print(f"#   best size match: {hits[0][0]} ({hits[0][1]}-{hits[0][2]} bp; {hits[0][3]})")
        if len(hits) > 1:
            others = ", ".join(h[0] for h in hits[1:4])
            print(f"#   other size-compatible: {others}")
    else:
        print(f"#   no cnvs.csv size match (novel/uncatalogued breakpoint)")
    print(f"#   VERDICT: DELETION ~{consensus}bp, {zyg}", file=sys.stderr)
    print(f"#     {args.name} ({args.genotype}) -> "
          f"{hits[0][0] if hits else 'uncatalogued'}", file=sys.stderr)

    zyg_short = ("HOM" if frac >= 0.85 else
                 "HET" if frac >= 0.30 else "UNCERTAIN")
    write_tsv(args.out, args.name, "yes", consensus,
              hits[0][0] if hits else "uncatalogued",
              zyg_short, f"{frac:.2f}", args.genotype)


if __name__ == "__main__":
    main()
