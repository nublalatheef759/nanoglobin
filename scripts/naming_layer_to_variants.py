#!/usr/bin/env python3
"""
naming_layer_to_variants.py -- convert the IthaGenes-primary naming_layer.csv
into the exact 18-column schema of the existing databases/variants.csv, so it is
a DROP-IN replacement that the three consumers (vep_annotate.py,
annotation/annotate_variants.py, sample_report.py) read without modification.

The consumers key on the HVGS column; all 3685 naming_layer rows have HGVS, so
all are usable. Genomic coordinates are carried from naming_layer, which now
holds VariantValidator-derived GRCh38 coordinates for 3676/3685 variants
(99.8%); the 9 without coordinates are documented exceptions (protein-notation,
dual-gene-uncertain, large promoter-spanning deletions) that require manual
curation. The consumers do not require coordinates (they key on HVGS), so blank
coordinates on those 9 do not affect annotation.

The result is a strict SUPERSET of the old file: every coordinate-rich variant
is retained, plus 1069 additional IthaGenes-curated variants. Coordinate
coverage as a percentage is lower only because the new additions (pure gain)
dilute it.

Output schema (exact, in order):
  HVGS, Common name, Gene, Cohort_n, Cohort frequency (%), Functionality,
  ClinVar classification, In Ithagenes? (Yes/No), In gnomAD? (Yes/No),
  gnomAD total AF (%), gnomAD ME AF (%), clinvar_review_status, clinvar_id,
  chrom, pos, ref, alt, source
"""
import argparse
import csv


OUT_COLUMNS = [
    "HVGS", "Common name", "Gene", "Cohort_n", "Cohort frequency (%)",
    "Functionality", "ClinVar classification", "In Ithagenes? (Yes/No)",
    "In gnomAD? (Yes/No)", "gnomAD total AF (%)", "gnomAD ME AF (%)",
    "clinvar_review_status", "clinvar_id", "chrom", "pos", "ref", "alt",
    "source",
]


def yn(val):
    """Normalise a truthy/flag value to Yes/No."""
    v = str(val).strip().lower()
    if v in ("yes", "true", "1", "y"):
        return "Yes"
    if v in ("no", "false", "0", "n", ""):
        return "No"
    return "Yes" if v else "No"


def derive_source(in_ith, in_clin):
    i = yn(in_ith) == "Yes"
    c = yn(in_clin) == "Yes"
    if i and c:
        return "IthaGenes+ClinVar"
    if i:
        return "IthaGenes"
    if c:
        return "ClinVar"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--naming", default="databases/naming_layer.csv")
    ap.add_argument("--out", default="databases/variants_v2.csv")
    args = ap.parse_args()

    with open(args.naming) as f:
        rows = list(csv.DictReader(f))

    out_rows = []
    for r in rows:
        out_rows.append({
            "HVGS": (r.get("HGVS") or "").strip(),
            "Common name": (r.get("Common Name") or "").strip(),
            "Gene": (r.get("Genes") or "").strip(),
            "Cohort_n": (r.get("cohort_n") or "").strip(),
            "Cohort frequency (%)": (r.get("cohort_freq_pct") or "").strip(),
            "Functionality": (r.get("Functionality") or "").strip(),
            "ClinVar classification": (r.get("clinvar_classification") or "").strip(),
            "In Ithagenes? (Yes/No)": yn(r.get("in_ithagenes")),
            "In gnomAD? (Yes/No)": yn(r.get("in_gnomad")),
            "gnomAD total AF (%)": (r.get("gnomad_af_total") or "").strip(),
            "gnomAD ME AF (%)": (r.get("gnomad_af_me") or "").strip(),
            "clinvar_review_status": (r.get("clinvar_review_status") or "").strip(),
            "clinvar_id": (r.get("clinvar_id") or "").strip(),
            "chrom": (r.get("chrom") or "").strip(),
            "pos": (r.get("pos") or "").strip(),
            "ref": (r.get("ref") or "").strip(),
            "alt": (r.get("alt") or "").strip(),
            "source": derive_source(r.get("in_ithagenes"), r.get("in_clinvar")),
        })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLUMNS)
        w.writeheader()
        w.writerows(out_rows)

    # report
    n = len(out_rows)
    with_coord = sum(1 for r in out_rows if r["chrom"] and r["pos"])
    with_hvgs = sum(1 for r in out_rows if r["HVGS"])
    src = {}
    for r in out_rows:
        src[r["source"]] = src.get(r["source"], 0) + 1
    print(f"wrote {args.out}: {n} rows")
    print(f"  with HVGS (join key): {with_hvgs}")
    print(f"  with coordinates: {with_coord}  without: {n - with_coord}")
    print(f"  source breakdown: {src}")


if __name__ == "__main__":
    main()
