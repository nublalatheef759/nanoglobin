#!/usr/bin/env python3
"""
detect_alpha_cnv.py -- copy-number detection for the HBA (alpha-globin) region
from Oxford Nanopore read depth.

This supersedes the earlier HBB-normalised / fixed-window approach, which two
findings during validation showed to be unreliable:

  1. HBB-normalisation is fragile: HBB coverage varies between samples
     (observed 20x-41x median across samples), so HBA/HBB ratios were distorted
     when a sample's HBB baseline was atypical. FIX: normalise HBA depth against
     a STABLE flanking region on the same chromosome (default chr16:1,000,000-
     1,100,000), which is unaffected by globin CNVs and tracks the sample's own
     sequencing depth.

  2. A fixed wide window dilutes localised deletions: a heterozygous -a4.2
     (~4.2 kb) measured over a 5 kb+ window that includes flanking normal
     sequence averages the 50% dip back up toward normal (observed 0.96 for a
     true het). FIX: scan the HBA region in small windows and report the
     DEEPEST contiguous dip, i.e. measure the deletion where it actually is.

Copy-number interpretation of the deepest-window ratio (test / normal-panel mean
at the same window), grounded in expected alpha-gene copy number (normal = 4):
    ~1.0        normal (4 copies)
    ~0.5        heterozygous single-gene / one-chromosome deletion (loses ~half
                the depth in the affected segment)
    ~0.0        homozygous deletion (both chromosomes affected -> no coverage)
    NOTE: absolute ratios are further depressed inside the HBA1/HBA2 segmental
    duplication because multi-mapping suppresses baseline depth even in normals;
    the het/hom DISTINCTION (residual vs zero coverage) is robust, the absolute
    value is region-dependent.

LIMITATION (documented, validated): alpha-TRIPLICATIONS (e.g. anti-3.7, aaa3.7)
are NOT detectable by read depth. The additional near-identical gene copy maps
to the same reference coordinates as the existing copies, so it produces no
depth increase (confirmed across mapq0, flanking-normalised and localised
analyses: an aaa3.7 carrier showed no elevated window anywhere). Triplication /
gain calling requires split-read / breakpoint-junction methods, not coverage.
"""
import argparse
import statistics as stats
import subprocess
import sys


def depth_windows(samtools, bam, chrom, start, end, win, mapq):
    """Return list of (win_start, mean_depth) over [start,end) in steps of win."""
    out = []
    for ws in range(start, end, win):
        we = min(ws + win, end)
        region = f"{chrom}:{ws}-{we}"
        proc = subprocess.run(
            [samtools, "depth", "-a", "-Q", str(mapq), "-r", region, bam],
            capture_output=True, text=True)
        if proc.returncode != 0:
            continue
        vals = [int(l.split("\t")[2]) for l in proc.stdout.splitlines() if l]
        if vals:
            out.append((ws, sum(vals) / len(vals)))
    return out


def region_mean(samtools, bam, chrom, start, end, mapq):
    region = f"{chrom}:{start}-{end}"
    proc = subprocess.run(
        [samtools, "depth", "-a", "-Q", str(mapq), "-r", region, bam],
        capture_output=True, text=True)
    vals = [int(l.split("\t")[2]) for l in proc.stdout.splitlines() if l]
    return (sum(vals) / len(vals)) if vals else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True, help="test (carrier) sorted BAM")
    ap.add_argument("--panel", nargs="+", required=True,
                    help="normal aa/aa sorted BAMs (the control panel)")
    ap.add_argument("--name", default="sample")
    ap.add_argument("--genotype", default="")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--hba", default="chr16:170000-178000",
                    help="HBA region to scan")
    ap.add_argument("--flank", default="chr16:1000000-1100000",
                    help="stable flanking region for self-normalisation")
    ap.add_argument("--win", type=int, default=500, help="scan window size (bp)")
    ap.add_argument("--mapq", type=int, default=20)
    args = ap.parse_args()

    hchrom, hspan = args.hba.split(":")
    hstart, hend = (int(x) for x in hspan.split("-"))
    fchrom, fspan = args.flank.split(":")
    fstart, fend = (int(x) for x in fspan.split("-"))

    # --- per-sample flanking baseline (own sequencing depth) ---
    test_flank = region_mean(args.samtools, args.bam, fchrom, fstart, fend, args.mapq)
    panel_flank = [region_mean(args.samtools, b, fchrom, fstart, fend, args.mapq)
                   for b in args.panel]

    # --- windowed HBA depth, flanking-normalised, for test and each panel member ---
    test_win = depth_windows(args.samtools, args.bam, hchrom, hstart, hend, args.win, args.mapq)
    panel_win = [depth_windows(args.samtools, b, hchrom, hstart, hend, args.win, args.mapq)
                 for b in args.panel]

    # panel mean (flanking-normalised) per window position
    panel_norm_by_pos = {}
    for pw, pf in zip(panel_win, panel_flank):
        for (ws, d) in pw:
            if pf > 0:
                panel_norm_by_pos.setdefault(ws, []).append(d / pf)
    panel_mean_by_pos = {ws: stats.mean(v) for ws, v in panel_norm_by_pos.items() if v}

    # test flanking-normalised, then ratio-to-panel per window
    print(f"# {args.name}" + (f"  (DRAGEN: {args.genotype})" if args.genotype else ""))
    print(f"# flanking baseline: test={test_flank:.2f}  "
          f"panel={','.join(f'{x:.1f}' for x in panel_flank)}")
    print("win_start\ttest_norm\tpanel_norm\tratio")
    ratios = []
    for (ws, d) in test_win:
        if test_flank <= 0 or ws not in panel_mean_by_pos:
            continue
        tnorm = d / test_flank
        pnorm = panel_mean_by_pos[ws]
        if pnorm <= 0:
            continue
        r = tnorm / pnorm
        ratios.append((ws, r, d))
        print(f"{ws}\t{tnorm:.3f}\t{pnorm:.3f}\t{r:.3f}")

    if not ratios:
        print("# no data", file=sys.stderr)
        return

    # --- deepest contiguous dip: find the min-ratio window and its neighbours ---
    rvals = [r for _, r, _ in ratios]
    min_ratio = min(rvals)
    min_ws = ratios[rvals.index(min_ratio)][0]
    # mean ratio over the contiguous run of low windows around the minimum
    # (windows within the affected block, ratio < 0.75)
    affected = [ws for ws, r, _ in ratios if r < 0.75]

    # --- verdict, grounded in copy-number expectation ---
    # (deepest window is the discriminator: ~0 hom, ~0.5 het, ~1 normal)
    if min_ratio <= 0.15:
        call = "HOMOZYGOUS deletion (deepest window ~0 -> both copies absent)"
    elif min_ratio <= 0.70:
        call = "HETEROZYGOUS / single-chromosome deletion (deepest window ~0.5)"
    else:
        call = "no deletion detected by depth"

    print("\n# ---- SUMMARY ----", file=sys.stderr)
    print(f"# {args.name} ({args.genotype})", file=sys.stderr)
    print(f"#   deepest window: chr16:{min_ws}  ratio={min_ratio:.3f}", file=sys.stderr)
    if affected:
        print(f"#   affected span: chr16:{min(affected)}-{max(affected)+args.win} "
              f"({len(affected)} windows < 0.75)", file=sys.stderr)
    print(f"#   VERDICT: {call}", file=sys.stderr)
    print(f"#   NOTE: triplications are not depth-detectable (see header); a "
          f"ratio near 1.0 does NOT exclude a gain.", file=sys.stderr)


if __name__ == "__main__":
    main()
