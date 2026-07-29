#!/usr/bin/env python3
"""
Panel-normalised CNV detection for the HBA region.

The HBA1/HBA2 segmental duplication depresses read-depth even in normal (aa/aa)
samples due to multi-mapping, so absolute coverage (or HBB-normalised coverage)
cannot cleanly reveal single-gene deletions. This script builds a per-bin
"expected normal" profile from a panel of aa/aa samples, then expresses each
test sample as (test / normal-panel-median) per bin.

Interpretation of the panel-normalised ratio, per bin:
  ~1.0        = normal copy number
  ~0.5-0.75   = heterozygous single-gene deletion (one of the region's copies lost)
  ~0.0-0.5    = homozygous / double-gene deletion (deeper loss)
  >1.25       = duplication / triplication (gain)

Because the baseline duplication effect is present in BOTH the panel and the
test sample, dividing cancels it out, leaving the true copy-number change.
"""
import argparse
import sys
import statistics as stats


def load_bins(path):
    """Return {(chrom,start): ratio_to_HBB} for HBA bins."""
    d = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6 or p[0] == "region" or p[0] != "HBA":
                continue
            chrom, start = p[1], int(p[2])
            try:
                ratio = float(p[5])   # ratio_to_HBB column
            except ValueError:
                continue
            d[(chrom, start)] = ratio
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", nargs="+", required=True,
                    help="coverage_bins.tsv of NORMAL aa/aa samples")
    ap.add_argument("--test", required=True,
                    help="coverage_bins.tsv of the test (carrier) sample")
    ap.add_argument("--name", default="test")
    ap.add_argument("--genotype", default="", help="DRAGEN genotype for reference")
    args = ap.parse_args()

    panel = [load_bins(p) for p in args.panel]
    test = load_bins(args.test)

    # per-bin normal median across the panel
    bins = sorted(test.keys())
    print(f"# Panel-normalised HBA coverage for {args.name}"
          + (f"  (DRAGEN: {args.genotype})" if args.genotype else ""))
    print(f"# panel = {len(panel)} normal samples")
    print("chrom\tstart\ttest_ratio\tnormal_median\tpanel_norm\tcall")

    norm_ratios = []
    for (chrom, start) in bins:
        normvals = [pn[(chrom, start)] for pn in panel if (chrom, start) in pn]
        if not normvals:
            continue
        normal_med = stats.median(normvals)
        if normal_med <= 0:
            continue
        pn_ratio = test[(chrom, start)] / normal_med
        norm_ratios.append(pn_ratio)
        # per-bin call
        if pn_ratio >= 1.25:
            call = "GAIN"
        elif pn_ratio <= 0.35:
            call = "HOM_DEL"
        elif pn_ratio <= 0.75:
            call = "HET_DEL"
        else:
            call = "."
        print(f"{chrom}\t{start}\t{test[(chrom,start)]:.3f}\t"
              f"{normal_med:.3f}\t{pn_ratio:.3f}\t{call}")

    # summary: find the largest CONTIGUOUS affected block (real CNVs are
    # contiguous runs, not scattered bins — scattered = noise)
    print("\n# ---- SUMMARY ----", file=sys.stderr)
    if norm_ratios:
        med = stats.median(norm_ratios)
        lo = min(norm_ratios)
        hi = max(norm_ratios)
        print(f"# {args.name}: median panel-norm ratio = {med:.3f} "
              f"(range {lo:.3f}-{hi:.3f})", file=sys.stderr)

        # classify each bin, then find longest contiguous run of each type
        def classify(r):
            if r >= 1.35:
                return "GAIN"
            if r <= 0.35:
                return "HOM_DEL"
            if r <= 0.75:
                return "HET_DEL"
            return "."
        labels = [classify(r) for r in norm_ratios]

        def longest_run(target):
            best = cur = 0
            for lb in labels:
                cur = cur + 1 if lb == target else 0
                best = max(best, cur)
            return best

        run_gain = longest_run("GAIN")
        run_hom = longest_run("HOM_DEL")
        run_het = longest_run("HET_DEL")
        n_gain = labels.count("GAIN")
        n_hom = labels.count("HOM_DEL")
        n_het = labels.count("HET_DEL")
        print(f"#   bins: GAIN={n_gain}(run {run_gain}) "
              f"HOM_DEL={n_hom}(run {run_hom}) "
              f"HET_DEL={n_het}(run {run_het})", file=sys.stderr)

        # verdict = largest CONTIGUOUS block wins (>=4 bins = >=800bp)
        MIN_RUN = 4
        runs = {"GAIN": run_gain, "HOM_DEL": run_hom, "HET_DEL": run_het}
        winner = max(runs, key=runs.get)
        if runs[winner] < MIN_RUN:
            v = "no clear contiguous copy-number change"
        elif winner == "GAIN":
            v = f"DUPLICATION/GAIN detected (contiguous run of {run_gain} bins)"
        elif winner == "HOM_DEL":
            v = f"HOMOZYGOUS/DOUBLE deletion detected (run of {run_hom} bins)"
        else:
            v = f"HETEROZYGOUS deletion detected (run of {run_het} bins)"
        print(f"#   VERDICT: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
