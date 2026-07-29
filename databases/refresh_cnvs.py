#!/usr/bin/env python3
"""
Rebuild cnvs.csv from the new IthaCNVs bulk CSV export (which post-dates, and
corrects, the previously HTML-scraped version).

The export has a messy 3-row header (grouped 5'/3' Breakpoint super-columns).
Real data starts at row 4 (0-indexed row 3). We map its chromosome-breakpoint
columns to the pipeline's cnvs.csv schema:

    name,chrom,start_min,start_max,end_min,end_max,size_min,size_max,
    exact,type,genes,functionality,locus,hgvs,ithaID,source

Breakpoint handling:
  - "Chromosome 5' BP from/to" -> start_min/start_max
  - "Chromosome 3' BP from/to" -> end_min/end_max
  - "No other BP defined" in the *to* field means an EXACT breakpoint:
    start_max := start_min (and likewise for end).
  - size = end_min - start_min (min) and end_max - start_max (max)

This is a STANDALONE artifact -> writes cnvs_refreshed.csv. It does NOT overwrite
your working cnvs.csv. Compare, then switch over deliberately if you choose.

Fields the export does NOT carry (functionality, hgvs, type=DEL/DUP) are filled
from your existing cnvs.csv by ithaID where possible, else left blank / inferred.
"""
import argparse
import csv
import sys
import pandas as pd
import numpy as np


NO_BP = "No other bp defined"


def clean(x):
    if x is None:
        return ""
    return str(x).strip()


def to_int(x):
    x = clean(x)
    if x == "" or x.lower() == NO_BP:
        return None
    try:
        return int(float(x.replace(",", "")))
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--export", default="ithaCNV_export_2026-07-28.csv")
    ap.add_argument("--old", default="cnvs.csv",
                    help="existing cnvs.csv to inherit functionality/hgvs/type from")
    ap.add_argument("--out", default="cnvs_refreshed.csv")
    ap.add_argument("--chrom", default="chr11:HBB,chr16:HBA",
                    help="not used for logic; chrom assigned per-gene below")
    args = ap.parse_args()

    # ---- read the export, skipping the 3 header rows; assign explicit columns ----
    # Column order (from the export header we inspected):
    # 0 IthaID, 1 Common name, 2 Genes,
    # 3 5'Type, 4 MLPA last-normal, 5 MLPA 1st-affected,
    # 6 3'... (MLPA last-affected), 7 MLPA 1st-normal, 8 3'Type,
    # 9 Locus5from,10 Locus5to,11 Locus3from,12 Locus3to,
    # 13 Chrom5from,14 Chrom5to,15 Chrom3from,16 Chrom3to
    cols = ["ithaID", "name", "genes",
            "bp5_type", "mlpa_last_normal", "mlpa_1st_affected",
            "mlpa_last_affected", "mlpa_1st_normal", "bp3_type",
            "locus5_from", "locus5_to", "locus3_from", "locus3_to",
            "chrom5_from", "chrom5_to", "chrom3_from", "chrom3_to"]

    rows = []
    with open(args.export, newline="", encoding="utf-8-sig") as fh:
        r = csv.reader(fh)
        all_rows = list(r)
    # data starts once IthaID column is numeric
    for raw in all_rows:
        if not raw:
            continue
        first = clean(raw[0])
        if first.isdigit():
            # pad/truncate to expected width
            raw = (raw + [""] * len(cols))[:len(cols)]
            rows.append(dict(zip(cols, [clean(c) for c in raw])))
    exp = pd.DataFrame(rows)
    print(f"Export data rows parsed: {len(exp)}", file=sys.stderr)

    # ---- inherit functionality/hgvs/type from old cnvs.csv by ithaID ----
    old = pd.read_csv(args.old, dtype=str, keep_default_na=False)
    old["ithaID"] = old["ithaID"].astype(str).str.strip()
    old_by_id = old.set_index("ithaID")
    def inherit(ithaid, field):
        try:
            return old_by_id.loc[str(ithaid), field]
        except (KeyError, TypeError):
            return ""

    # ---- build the refreshed rows in pipeline schema ----
    out_rows = []
    for _, e in exp.iterrows():
        ithaid = e["ithaID"]
        genes = e["genes"]

        # chromosome: HBB cluster = chr11, HBA cluster = chr16
        g = genes.upper()
        if any(k in g for k in ["HBA", "HBZ", "HS40", "HBM", "HBQ"]):
            chrom = "chr16"
        else:
            chrom = "chr11"   # HBB/HBD/HBG/HBE cluster

        start_from = to_int(e["chrom5_from"])
        start_to = to_int(e["chrom5_to"])
        end_from = to_int(e["chrom3_from"])
        end_to = to_int(e["chrom3_to"])

        # Breakpoints may populate either the 'from' or 'to' column (the other
        # reads "No other BP defined"). Collect whatever numeric values exist
        # and take min/max, so exact breakpoints (single value in either column)
        # and true ranges (both columns) both work.
        start_vals = [v for v in (start_from, start_to) if v is not None]
        end_vals = [v for v in (end_from, end_to) if v is not None]
        start_min = min(start_vals) if start_vals else None
        start_max = max(start_vals) if start_vals else None
        end_min = min(end_vals) if end_vals else None
        end_max = max(end_vals) if end_vals else None

        # exactness from the Type fields
        exact = "yes" if (clean(e["bp5_type"]).lower() == "exact"
                          and clean(e["bp3_type"]).lower() == "exact") else "no"

        # sizes
        def size(a, b):
            if a is None or b is None:
                return ""
            return abs(b - a)
        size_min = size(start_min, end_min)
        size_max = size(start_max, end_max)

        # inherited fields
        functionality = inherit(ithaid, "functionality")
        hgvs = inherit(ithaid, "hgvs")
        old_type = inherit(ithaid, "type")
        locus = inherit(ithaid, "locus")
        # infer type if not inherited: default DEL (IthaCNVs is deletions/dups)
        vtype = old_type if old_type else "DEL"

        out_rows.append({
            "name": e["name"],
            "chrom": chrom,
            "start_min": start_min if start_min is not None else "",
            "start_max": start_max if start_max is not None else "",
            "end_min": end_min if end_min is not None else "",
            "end_max": end_max if end_max is not None else "",
            "size_min": size_min,
            "size_max": size_max,
            "exact": exact,
            "type": vtype,
            "genes": genes,
            "functionality": functionality,
            "locus": locus,
            "hgvs": hgvs,
            "ithaID": ithaid,
            "source": f"IthaCNVs ithaID={ithaid} (GRCh38.p13, export 2026-07-28)",
        })

    out = pd.DataFrame(out_rows, columns=[
        "name", "chrom", "start_min", "start_max", "end_min", "end_max",
        "size_min", "size_max", "exact", "type", "genes", "functionality",
        "locus", "hgvs", "ithaID", "source"])
    out.to_csv(args.out, index=False)

    # ---- report: what changed vs old ----
    old_ids = set(old["ithaID"])
    new_ids = set(out["ithaID"])
    added = new_ids - old_ids
    removed = old_ids - new_ids
    print("=" * 55, file=sys.stderr)
    print(f"cnvs_refreshed.csv written: {len(out)} CNVs", file=sys.stderr)
    print(f"  in old cnvs.csv:          {len(old)}", file=sys.stderr)
    print(f"  NEW (added in export):    {len(added)}  {sorted(added)[:20]}", file=sys.stderr)
    print(f"  in old but NOT in export: {len(removed)}  {sorted(removed)[:20]}", file=sys.stderr)
    print("=" * 55, file=sys.stderr)

    # coordinate diffs on shared IDs (did any breakpoints change? = the bug fixes)
    print("Coordinate changes on shared variants (start_min/end_min):", file=sys.stderr)
    changed = 0
    for iid in sorted(new_ids & old_ids, key=lambda x: (len(x), x)):
        try:
            o = old_by_id.loc[iid]
            n = out[out["ithaID"] == iid].iloc[0]
            if str(o["start_min"]) != str(n["start_min"]) or str(o["end_min"]) != str(n["end_min"]):
                print(f"  ithaID {iid} ({n['name'][:30]}): "
                      f"start {o['start_min']}->{n['start_min']}, "
                      f"end {o['end_min']}->{n['end_min']}", file=sys.stderr)
                changed += 1
        except Exception:
            pass
    if changed == 0:
        print("  (none — coordinates identical for all shared variants)", file=sys.stderr)
    else:
        print(f"  {changed} variants with changed coordinates", file=sys.stderr)


if __name__ == "__main__":
    main()
