#!/usr/bin/env python3
"""
concordance.py -- add cross-caller concordance to the comprehensive report.

For every SV row (Sniffles or CuteSV), check whether the OTHER caller made an
overlapping call in the same sample (reciprocal overlap >= 0.5 AND size within
+/-20%). Adds two columns:
  Concordance  : "both" if the other caller agrees, else "single"
  Concordant_with : the matching call's position from the other caller (or "")

Two independent callers agreeing is stronger evidence than either alone. This is
a confidence annotation, not a filter -- single-caller calls are kept (dropping
them would lose real variants the other caller happened to miss), just marked.

Position/size jitter between callers is expected (e.g. Sniffles 173707/-3812 vs
CuteSV 173708/-3811 for the same -a3.7), so matching is by overlap, not identity.

Usage:
  python concordance.py results/comprehensive_report.csv results/comprehensive_report.csv
"""
import csv, re, sys

RECIP, SIZE_TOL = 0.5, 0.20
SVLEN_RE = re.compile(r"SV_[A-Z]+_(-?\d+)bp")


def span(row):
    """(chrom, start, end) for an SV report row, size from the Consequence field."""
    m = SVLEN_RE.search(row.get("Consequence", ""))
    if not m:
        return None
    start = int(row["Position"])
    size = abs(int(m.group(1)))
    return row["Chromosome"], start, start + size


def overlaps(a, b):
    if a is None or b is None or a[0] != b[0]:
        return False
    (_, s1, e1), (_, s2, e2) = a, b
    l1, l2 = e1 - s1, e2 - s2
    if l1 <= 0 or l2 <= 0:
        return False
    ov = max(0, min(e1, e2) - max(s1, s2))
    if min(ov / l1, ov / l2) < RECIP:
        return False
    return (1 - SIZE_TOL) <= (l2 / l1) <= (1 + SIZE_TOL)


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    with open(inp, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys()) if rows else []

    for c in ("Concordance", "Concordant_with"):
        if c not in fields:
            fields.append(c)

    # index SV calls by sample+tool
    sv = {}  # (sample, tool) -> list of (row_index, span)
    for i, r in enumerate(rows):
        if r["Tool"] in ("Sniffles", "CuteSV"):
            sp = span(r)
            if sp:
                sv.setdefault((r["Sample"], r["Tool"]), []).append((i, sp))

    for i, r in enumerate(rows):
        if r["Tool"] not in ("Sniffles", "CuteSV"):
            r["Concordance"] = ""       # SNVs: not applicable
            r["Concordant_with"] = ""
            continue
        other = "CuteSV" if r["Tool"] == "Sniffles" else "Sniffles"
        mine = span(r)
        match = ""
        for j, sp in sv.get((r["Sample"], other), []):
            if overlaps(mine, sp):
                match = "%s:%d" % (rows[j]["Chromosome"], int(rows[j]["Position"]))
                break
        r["Concordance"] = "both" if match else "single"
        r["Concordant_with"] = match

    with open(outp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # summary
    both = sum(1 for r in rows if r.get("Concordance") == "both")
    single = sum(1 for r in rows if r.get("Concordance") == "single")
    sys.stderr.write("SV calls: %d concordant (both callers), %d single-caller\n" % (both, single))


if __name__ == "__main__":
    main()
    
    