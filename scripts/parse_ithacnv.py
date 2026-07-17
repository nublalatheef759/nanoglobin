#!/usr/bin/env python3
"""
parse_ithacnv.py -- turn the IthaCNVs HTML table into databases/cnvs.csv

Source: https://www.ithanet.eu/db/ithacnv  (coordinates are GRCh38.p13)
  wget -O ithacnv.html "https://www.ithanet.eu/db/ithacnv"
  python scripts/parse_ithacnv.py ithacnv.html databases/cnvs.csv

IthaCNVs reports each breakpoint as a pair (from, to):
  - "Exact" breakpoints put the position in one field and "No other BP defined"
    in the other, i.e. start = 5' BP from, end = 3' BP to.
  - "Range" breakpoints populate both, because the junction was never sequenced
    and is only localised to an interval.
  - "Unknown" entries say "Information unclear" and are kept but left blank.

The output preserves that uncertainty explicitly:
  start_min, start_max, end_min, end_max
For an exact breakpoint start_min == start_max. Downstream matching should use
the interval, not pretend a range is a point.
"""

import re
import sys

import pandas as pd

NON_NUMERIC = (
    "no other bp defined", "information unclear", "out of locus", "unknown",
    "bp extends beyond the mlpa probes", "bp starts from the telomere", "nan", "",
)


def flatten(df):
    """Collapse the 3-level MultiIndex header to a single meaningful name per column."""
    names = []
    for col in df.columns:
        parts = [str(p) for p in col if not str(p).startswith("Unnamed")]
        # de-duplicate repeated level labels e.g. ('IthaID','IthaID','IthaID')
        seen, uniq = set(), []
        for p in parts:
            if p not in seen:
                seen.add(p)
                uniq.append(p)
        names.append(" | ".join(uniq) if uniq else "blank")
    df = df.copy()
    df.columns = names
    return df


def find_col(df, *musts):
    """First column whose name contains all the given substrings (case-insensitive)."""
    for c in df.columns:
        lc = c.lower()
        if all(m.lower() in lc for m in musts):
            return c
    return None


def num(v):
    """Return int, or None for any of IthaCNVs' many non-numeric placeholders."""
    if v is None:
        return None
    s = str(v).strip()
    if s.lower() in NON_NUMERIC:
        return None
    s = s.replace(",", "").replace(" ", "")
    return int(s) if re.fullmatch(r"\d+", s) else None


def sv_type(hgvs):
    """Infer SV type from the HGVS string."""
    h = (hgvs or "").lower()
    if "delins" in h:
        return "DELINS"
    if "dup" in h and "del" in h:
        return "COMPLEX"
    if "dup" in h:
        return "DUP"
    if "del" in h:
        return "DEL"
    return ""


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: parse_ithacnv.py <ithacnv.html> <out.csv>")
    html, out = sys.argv[1], sys.argv[2]

    tables = pd.read_html(html)
    df = max(tables, key=lambda t: t.shape[0])       # the 311-row one
    if isinstance(df.columns, pd.MultiIndex):
        df = flatten(df)

    cols = {
        "itha":      find_col(df, "IthaID"),
        "name":      find_col(df, "Common name"),
        "genes":     find_col(df, "Genes"),
        "hgvs":      find_col(df, "HGVS name"),
        "func":      find_col(df, "Functionality"),
        "chrom":     "Chromosome" if "Chromosome" in df.columns else None,
        "s_from":    find_col(df, "Chromosome 5' BP from"),
        "s_to":      find_col(df, "Chromosome 5' BP to"),
        "e_from":    find_col(df, "Chromosome 3' BP from"),
        "e_to":      find_col(df, "Chromosome 3' BP to"),
        "locus":     "Locus" if "Locus" in df.columns else None,
    }
    missing = [k for k, v in cols.items() if v is None and k not in ("chrom", "locus")]
    if missing:
        sys.stderr.write("Columns found:\n")
        for c in df.columns:
            sys.stderr.write("  %s\n" % c)
        sys.exit("ERROR: could not locate columns: %s" % missing)

    rows, skipped = [], 0
    for _, r in df.iterrows():
        itha = num(r.get(cols["itha"]))
        if itha is None:
            continue

        s_from, s_to = num(r.get(cols["s_from"])), num(r.get(cols["s_to"]))
        e_from, e_to = num(r.get(cols["e_from"])), num(r.get(cols["e_to"]))

        # exact: one side of the pair is a placeholder -> collapse to a point
        start_min = s_from
        start_max = s_to if s_to is not None else s_from
        end_max = e_to
        end_min = e_from if e_from is not None else e_to

        if start_min is None or end_max is None:
            skipped += 1
            continue
        if end_max <= start_min:
            sys.stderr.write("WARNING ithaID=%s %s: start %s > end %s, skipping\n"
                             % (itha, r.get(cols["name"]), start_min, end_max))
            skipped += 1
            continue
        if start_max is not None and start_max < start_min:
            start_min, start_max = start_max, start_min
        if end_min is not None and end_min > end_max:
            end_min, end_max = end_max, end_min

        chrom_raw = str(r.get(cols["chrom"])).strip() if cols["chrom"] else ""
        chrom = "chr" + chrom_raw if re.fullmatch(r"\d+", chrom_raw) else ""
        hgvs = str(r.get(cols["hgvs"]) or "").strip()

        rows.append({
            "name": str(r.get(cols["name"]) or "").strip(),
            "chrom": chrom,
            "start_min": start_min,
            "start_max": start_max,
            "end_min": end_min,
            "end_max": end_max,
            "size_min": max(0, end_min - start_max) if end_min else "",
            "size_max": end_max - start_min,
            "exact": "yes" if (start_min == start_max and end_min == end_max) else "no",
            "type": sv_type(hgvs),
            "genes": str(r.get(cols["genes"]) or "").strip(),
            "functionality": str(r.get(cols["func"]) or "").strip(),
            "locus": str(r.get(cols["locus"]) or "").strip() if cols["locus"] else "",
            "hgvs": hgvs,
            "ithaID": itha,
            "source": "IthaCNVs ithaID=%d (GRCh38.p13)" % itha,
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out, index=False)

    sys.stderr.write("parsed %d rows, skipped %d with unusable breakpoints\n"
                     % (len(out_df), skipped))
    sys.stderr.write("wrote %s\n\n" % out)
    sys.stderr.write("by chromosome:\n%s\n\n" % out_df["chrom"].value_counts().to_string())
    sys.stderr.write("by type:\n%s\n\n" % out_df["type"].value_counts().to_string())
    sys.stderr.write("exact vs range:\n%s\n" % out_df["exact"].value_counts().to_string())


if __name__ == "__main__":
    main()
    