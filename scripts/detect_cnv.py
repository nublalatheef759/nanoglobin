#!/usr/bin/env python3
"""
detect_cnv.py -- call copy-number GAINS and LOSSES from the binned coverage
profile (the *_coverage_bins.tsv produced by coverage_profile.py).

The whole-window mean masks localised events through flanking dilution (a 3.8kb
gain in an 8kb window barely moves the mean). This works on the binned profile
instead: it finds a run of consecutive bins that departs from the local flanking
baseline, which recovers both the event and its boundaries.

  loss  : run of bins <= 0.5 * flanking baseline   (deletion)
  gain  : run of bins >= 1.5 * flanking baseline   (duplication/triplication)

Baseline = median of the bins OUTSIDE the called run (the flanking normal depth),
so it is robust to the HBA/HBB ratio not being exactly 1.0.
"""
import csv, sys, statistics

MIN_RUN = 3          # consecutive bins required
GAIN_FACTOR = 1.5    # >= 1.5x baseline -> gain
LOSS_FACTOR = 0.5    # <= 0.5x baseline -> loss


def load_bins(path, region="HBA"):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["region"] == region:
                rows.append((int(r["bin_start"]), float(r["ratio_to_HBB"])))
    rows.sort()
    return rows


# The ratio is already HBB-normalised, so 1.0 IS the diploid baseline. Do NOT
# use the sample's own median: a whole-region event (e.g. flat 0.5 for full-alpha)
# moves the median and would hide itself. Fixed baseline of 1.0 avoids that.
BASELINE = 1.0


def find_runs(rows, baseline=BASELINE):
    calls = []
    def classify(ratio):
        if baseline <= 0 or ratio != ratio:
            return None
        f = ratio / baseline
        if f >= GAIN_FACTOR: return "gain"
        if f <= LOSS_FACTOR: return "loss"
        return None
    i, n = 0, len(rows)
    while i < n:
        kind = classify(rows[i][1])
        if kind is None:
            i += 1; continue
        j = i
        while j < n and classify(rows[j][1]) == kind:
            j += 1
        if j - i >= MIN_RUN:
            seg = rows[i:j]
            start = seg[0][0]
            end = seg[-1][0]
            med = statistics.median([r for _, r in seg])
            copies = round(2 * med / baseline) if baseline else None
            calls.append((kind, start, end, med / baseline, copies))
        i = j
    return calls


def call_file(path):
    rows = load_bins(path)
    if not rows:
        return []
    return find_runs(rows, BASELINE), BASELINE


def write_cnv_calls(out_path, samples):
    """Write one row per CNV call to a CSV the sample report can consume.

    Gains are matched against the IthaCNVs catalogue by reciprocal overlap
    (reusing identify_sv) so a detected bump can be named, e.g. ααα(anti-3.7).
    Copy number for gains is noisy (a duplication maps ambiguously), so we report
    the measured fold-change, not an inferred copy count."""
    import glob, os, csv as _csv, sys as _sys
    _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from identify_sv import identify_sv
    except Exception:
        identify_sv = None
    if not samples:
        _sys.exit("ERROR: no samples declared. "
                  "usage: detect_cnv.py <out.csv> <sample> [sample ...]")
    paths = []
    for s in samples:
        hits = sorted(glob.glob("variants/%s/*.coverage_bins.tsv" % s))
        if not hits:
            _sys.exit("ERROR: no coverage bins found for declared sample %s" % s)
        paths.extend(hits)
    with open(out_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = _csv.writer(fh)
        w.writerow(["Sample", "Type", "Region", "Fold", "Name", "Confidence"])
        for p in paths:
            sample = os.path.basename(os.path.dirname(p))
            res = call_file(p)
            if not res:
                continue
            calls, base = res
            for kind, start, end, fold, copies in calls:
                region = "chr16:%d-%d" % (start, end)
                name = ""
                if kind == "gain" and identify_sv is not None:
                    # match the bump against catalogued DUP entries
                    hit = identify_sv("chr16", start, "DUP", end - start)
                    if hit and not hit.startswith("unknown"):
                        name = hit
                conf = "high" if kind == "loss" else "moderate"
                w.writerow([sample, kind, region, "%.2f" % fold, name, conf])


if __name__ == "__main__":
    import glob, os
    # test across all samples: gains should fire only on triple, losses on deletions
    paths = sorted(glob.glob("variants/*/*.coverage_bins.tsv"))
    print("%-16s %-6s %s" % ("sample", "base", "calls"))
    for p in paths:
        sample = os.path.basename(os.path.dirname(p))
        res = call_file(p)
        if not res:
            print("%-16s   --   (no HBA bins)" % sample); continue
        calls, base = res
        if not calls:
            print("%-16s %.2f   normal (no CNV)" % (sample, base))
        else:
            desc = "; ".join("%s chr16:%d-%d (%.2fx, ~%s copies)" % (k, s, e, f, c)
                             for k, s, e, f, c in calls)
            print("%-16s %.2f   %s" % (sample, base, desc))

    import sys as _sys
    if len(_sys.argv) > 1:
        if len(_sys.argv) < 3:
            _sys.exit("usage: detect_cnv.py <out.csv> <sample> [sample ...]")
        write_cnv_calls(_sys.argv[1], _sys.argv[2:])
        print("\nwrote CNV calls -> %s" % _sys.argv[1])
        
        
            