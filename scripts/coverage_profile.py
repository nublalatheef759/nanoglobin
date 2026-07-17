#!/usr/bin/env python3
"""
coverage_profile.py -- binned read depth across target regions, normalised to a
reference region on a different chromosome.

Why normalise against another chromosome:
  Comparing a region's depth against itself cannot detect a deletion that removes
  the whole region: every position drops by the same factor, the internal ratios
  are unchanged, and the region looks normal. This is the documented failure mode
  of within-sample normalisation for full alpha-cluster deletions (--/aa), which
  get reported as aa/aa. Normalising HBA against HBB -- a locus on chr11 that the
  deletion does not touch -- gives a ratio that halves.

Why bin rather than take a single mean:
  A single mean over an 8kb window catches --/aa (whole window drops) but cannot
  localise a 3.8kb deletion inside it, and cannot distinguish -a3.7 from -a4.2,
  which remove different genes. Binned depth turns a deletion into a shape:
  normal, drop, normal -- with the drop where the deleted segment is.

samtools depth is called with -a so that zero-coverage positions are emitted.
Without -a those positions are silently omitted and a homozygous deletion is
averaged out of existence.

Usage:
  python coverage_profile.py \
      --bam sorted_reads/S.bam \
      --target HBA=chr16:170000-178000 \
      --target HBB=chr11:5225000-5228000 \
      --reference HBB \
      --bin-size 200 \
      --summary variants/S/S.coverage.tsv \
      --bins variants/S/S.coverage_bins.tsv
"""

import argparse
import statistics
import subprocess
import sys


def parse_region(spec):
    """'HBA=chr16:170000-178000' -> ('HBA', 'chr16', 170000, 178000)"""
    name, _, loc = spec.partition("=")
    if not name or not loc:
        raise ValueError("bad --target %r, expected NAME=chrom:start-end" % spec)
    chrom, _, span = loc.partition(":")
    start, _, end = span.partition("-")
    return name, chrom, int(start), int(end)


def samtools_depth(bam, chrom, start, end, samtools="samtools"):
    """Return a list of depths, one per base from start to end inclusive.

    -a forces output of zero-coverage positions. Without it, deleted regions
    produce no rows at all and vanish from any downstream mean.
    """
    region = "%s:%d-%d" % (chrom, start, end)
    proc = subprocess.run(
        [samtools, "depth", "-a", "-r", region, bam],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        sys.exit("samtools depth failed for %s:\n%s" % (region, proc.stderr))

    by_pos = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            by_pos[int(f[1])] = int(f[2])

    # samtools omits positions past the end of a contig; treat those as zero
    return [by_pos.get(p, 0) for p in range(start, end + 1)]


def bin_depths(depths, start, bin_size):
    """Yield (bin_start, bin_end, mean_depth) over the depth array."""
    for i in range(0, len(depths), bin_size):
        chunk = depths[i:i + bin_size]
        if not chunk:
            continue
        yield start + i, start + i + len(chunk) - 1, sum(chunk) / len(chunk)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--target", action="append", required=True,
                   help="NAME=chrom:start-end, repeatable")
    p.add_argument("--reference", required=True,
                   help="name of the target to normalise against")
    p.add_argument("--bin-size", type=int, default=200)
    p.add_argument("--summary", required=True)
    p.add_argument("--bins", required=True)
    p.add_argument("--samtools", default="samtools")
    a = p.parse_args()

    regions = [parse_region(t) for t in a.target]
    names = [r[0] for r in regions]
    if a.reference not in names:
        sys.exit("--reference %s is not one of the targets: %s" % (a.reference, names))

    depths = {}
    for name, chrom, start, end in regions:
        depths[name] = (chrom, start, samtools_depth(a.bam, chrom, start, end, a.samtools))

    # Median, not mean, for the denominator: a median is unmoved by a deletion
    # affecting part of the reference region, and by amplicon edge effects.
    ref_chrom, ref_start, ref_depths = depths[a.reference]
    ref_median = statistics.median(ref_depths) if ref_depths else 0.0
    if ref_median == 0:
        sys.stderr.write(
            "WARNING: reference region %s has median depth 0. Ratios will be "
            "undefined. Check the BAM covers %s.\n" % (a.reference, ref_chrom)
        )

    # ---- per-bin profile -------------------------------------------------
    with open(a.bins, "w") as fh:
        fh.write("region\tchrom\tbin_start\tbin_end\tmean_depth\tratio_to_%s\n" % a.reference)
        for name, chrom, start, end in regions:
            _, _, d = depths[name]
            for b_start, b_end, mean in bin_depths(d, start, a.bin_size):
                ratio = (mean / ref_median) if ref_median else float("nan")
                fh.write("%s\t%s\t%d\t%d\t%.3f\t%.4f\n"
                         % (name, chrom, b_start, b_end, mean, ratio))

    # ---- summary ---------------------------------------------------------
    # Keeps the original two <NAME>_mean_depth lines so comprehensive_report.py
    # continues to parse it, and appends the ratios.
    with open(a.summary, "w") as fh:
        for name, chrom, start, end in regions:
            _, _, d = depths[name]
            mean = sum(d) / len(d) if d else 0.0
            fh.write("%s_mean_depth\t%.3f\n" % (name, mean))
        for name, chrom, start, end in regions:
            if name == a.reference:
                continue
            _, _, d = depths[name]
            med = statistics.median(d) if d else 0.0
            ratio = (med / ref_median) if ref_median else float("nan")
            fh.write("%s_%s_ratio\t%.4f\n" % (name, a.reference, ratio))
        fh.write("%s_median_depth\t%.3f\n" % (a.reference, ref_median))

    sys.stderr.write("reference %s median depth: %.1f\n" % (a.reference, ref_median))
    for name, chrom, start, end in regions:
        _, _, d = depths[name]
        med = statistics.median(d) if d else 0.0
        sys.stderr.write("  %-6s median %9.1f   ratio %.3f\n"
                         % (name, med, (med / ref_median) if ref_median else float("nan")))


if __name__ == "__main__":
    main()
    