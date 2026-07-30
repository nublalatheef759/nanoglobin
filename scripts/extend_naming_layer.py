#!/usr/bin/env python3
"""
extend_naming_layer.py -- add HBD/HBG1/HBG2/HBE1 variants from the full IthaGenes
export into naming_layer.csv, extending the catalogue beyond HBA/HBB.

The consumers (sample_report.py, annotate_variants.py) key on the HGVS column,
so these rows are usable for naming without genomic coordinates (coord columns
left blank, exactly as 746 existing rows already are). coord_status records why.

Idempotent: skips genes already present; refuses to add duplicate HGVS.

Usage:
  python3 extend_naming_layer.py \
      --full databases/IthaGenes_full.csv \
      --naming databases/naming_layer.csv \
      --genes HBD HBG1 HBG2 HBE1
Writes naming_layer.csv in place (back it up first).
"""
import argparse, csv, sys, re

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", default="databases/IthaGenes_full.csv")
    ap.add_argument("--naming", default="databases/naming_layer.csv")
    ap.add_argument("--genes", nargs="+", default=["HBD", "HBG1", "HBG2", "HBE1"])
    args = ap.parse_args()

    # --- read the target schema (naming_layer header, exact order) ---
    with open(args.naming, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
    n_cols = len(header)
    # index the header so we place values by NAME not position
    hidx = {h: i for i, h in enumerate(header)}

    # --- read existing HGVS to avoid duplicates ---
    with open(args.naming, newline="") as f:
        existing = set()
        r = csv.DictReader(f)
        hgvs_col = header[0]  # first column is HGVS
        for row in r:
            existing.add((row.get(hgvs_col) or "").strip())

    # --- pull matching rows from the full IthaGenes export ---
    gene_pat = re.compile(r"^(" + "|".join(re.escape(g) for g in args.genes) + r"):c", re.I)
    added = []
    with open(args.full, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            hgvs = (row.get("HGVS Name") or "").strip()
            if not gene_pat.match(hgvs):
                continue
            if hgvs in existing:
                continue
            # build a full-width row keyed by the naming_layer header
            out = [""] * n_cols
            def put(col, val):
                if col in hidx:
                    out[hidx[col]] = val
            put(header[0], hgvs)                                  # HGVS
            put("Common Name", (row.get("Common Name") or "").strip())
            put("Hb Name", (row.get("Hb Name") or "").strip())
            put("Genes", (row.get("Genes") or "").strip())
            put("Functionality", (row.get("Functionality") or "").strip())
            put("Phenotype", (row.get("Phenotype") or "").strip())
            put("Locus", (row.get("Locus") or "").strip())
            put("Position", (row.get("Position") or "").strip())
            put("IthaID", (row.get("IthaID") or "").strip())
            # IthaGenes-only provenance; no ClinVar match asserted
            put("in_ithagenes", "Yes")
            put("in_clinvar", "No")
            put("coord_status", "not_derived_hgvs_keyed_extension")
            # clinvar_*, gnomad_*, cohort_*, chrom/pos/ref/alt left blank
            added.append(out)
            existing.add(hgvs)

    if not added:
        print("No new variants to add (already present or none matched).")
        return

    # --- append ---
    with open(args.naming, "a", newline="") as f:
        w = csv.writer(f)
        w.writerows(added)

    # --- report ---
    by_gene = {}
    for out in added:
        h = out[0]
        g = h.split(":")[0]
        by_gene[g] = by_gene.get(g, 0) + 1
    print(f"Added {len(added)} variants to {args.naming}:")
    for g, n in sorted(by_gene.items()):
        print(f"  {g}: {n}")
    print("Coordinates left blank (consumers key on HGVS). "
          "coord_status=not_derived_hgvs_keyed_extension")

if __name__ == "__main__":
    main()
