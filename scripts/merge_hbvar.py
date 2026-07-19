#!/usr/bin/env python3
"""
merge_hbvar.py -- add HbVar common names to the catalogue.

HbVar's tab-separated export gives common names (Hb Antananarivo, IVS I-110)
keyed by HGVS. We use it to fill the Common name field for catalogue entries
(mostly ClinVar-derived) that have a classification but no name.

HbVar quirks handled:
  - every field is wrapped in single quotes: 'Hb Sawara'
  - HGVS carries a paralogue hedge: 'HBA2:c.20A>C (or HBA1)' -> HBA2:c.20A>C
  - many hgvs_name values are not GENE:c. notation (NG_..:g., rsIDs, prose);
    those cannot join and are skipped.

Join key is GENE:c.notation, matching the catalogue's HVGS column.
Existing catalogue names are NEVER overwritten -- HbVar only fills blanks.

Usage:
  python merge_hbvar.py --hbvar hbvar.tsv --catalogue databases/variants.csv \
      --out databases/variants.csv
"""
import argparse, csv, re, sys

def unq(s):
    s = (s or "").strip()
    if len(s) >= 2 and s[0] == "'" and s[-1] == "'":
        s = s[1:-1]
    return s.strip()

CLEAN_HGVS = re.compile(r"^(HB[ABDGEZ][12]?):(c\.[^\s]+)")

def hbvar_names(path):
    """Yield (GENE:c.notation, common_name) for parseable rows."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        cols = {h.lstrip("#").strip(): i for i, h in enumerate(header)}
        hi, ni = cols.get("hgvs_name"), cols.get("hb_name")
        si = cols.get("syll_name")
        if hi is None or ni is None:
            sys.exit("could not find hgvs_name / hb_name columns")
        for row in rd:
            if len(row) <= max(hi, ni):
                continue
            hgvs = unq(row[hi]).replace(" (or HBA1)", "").replace(" (or HBA2)", "")
            m = CLEAN_HGVS.match(hgvs)
            if not m:
                continue
            key = "%s:%s" % (m.group(1), m.group(2))
            name = unq(row[ni])
            # prefer a real Hb name; fall back to syllabic (IVS/codon) name
            if not name or name.upper() in ("NULL", ""):
                name = unq(row[si]) if si is not None and si < len(row) else ""
            if name and name.upper() != "NULL":
                yield key, name

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hbvar", required=True)
    ap.add_argument("--catalogue", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    names = {}
    for key, name in hbvar_names(a.hbvar):
        names.setdefault(key, name)     # first wins
    sys.stderr.write("HbVar parseable GENE:c. names: %d\n" % len(names))

    with open(a.catalogue, newline="", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh)
        fields = rd.fieldnames
        rows = list(rd)

    filled = 0
    for r in rows:
        if not (r.get("Common name") or "").strip():
            nm = names.get((r.get("HVGS") or "").strip())
            if nm:
                r["Common name"] = nm
                filled += 1

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    have_name = sum(1 for r in rows if (r.get("Common name") or "").strip())
    sys.stderr.write("catalogue rows:            %d\n" % len(rows))
    sys.stderr.write("names filled from HbVar:   %d\n" % filled)
    sys.stderr.write("total with common name:    %d\n" % have_name)
    sys.stderr.write("wrote %s\n" % a.out)

if __name__ == "__main__":
    main()
    