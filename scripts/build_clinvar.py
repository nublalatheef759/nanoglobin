#!/usr/bin/env python3
"""
build_clinvar.py -- merge ClinVar globin-gene variants into the curated catalogue,
producing databases/variants.csv: a universal naming/classification layer that
works for any lab's data, not just the source cohort.

The curated catalogue (variant_lookup.csv) supplies common names and cohort data.
ClinVar supplies breadth and pathogenicity for everything else. Where both know a
variant, the curated common name is kept and ClinVar's classification is added.

Match key is GENE:c.notation (transcript-independent), so it is stable regardless
of which transcript VEP happens to report.

Input:  ClinVar variant_summary.txt.gz
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

Usage:
  python build_clinvar.py --clinvar variant_summary.txt.gz \
      --curated databases/variant_lookup.csv --out databases/variants.csv
"""
import argparse, csv, gzip, re, sys

GLOBIN = {"HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "HBE1", "HBZ"}
NAME_RE = re.compile(r"\(([A-Z0-9]+)\):(c\.[^ ]+)")   # ...(HBB):c.20A>T...


def clinvar_rows(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("GeneSymbol") not in GLOBIN:
                continue
            if row.get("Assembly") != "GRCh38":
                continue
            m = NAME_RE.search(row.get("Name", ""))
            if not m:
                continue
            gene, c = m.group(1), m.group(2)
            if gene not in GLOBIN:
                continue
            yield "%s:%s" % (gene, c), gene, row


def load_curated(path):
    by_key = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh)
        rd.fieldnames = [(f or "").strip() for f in rd.fieldnames]
        for row in rd:
            row = {(k or "").strip(): (v or "").strip() for k, v in row.items() if k}
            h = row.get("HVGS", "").strip()
            if h:
                by_key[h] = row
    return by_key


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinvar", required=True)
    ap.add_argument("--curated", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    curated = load_curated(a.curated)
    clin = {}
    for key, gene, row in clinvar_rows(a.clinvar):
        prev = clin.get(key)
        if prev and int(row.get("NumberSubmitters") or 0) <= int(prev[1]):
            continue
        clin[key] = (row, int(row.get("NumberSubmitters") or 0), gene)

    sys.stderr.write("curated variants:            %d\n" % len(curated))
    sys.stderr.write("ClinVar globin (GRCh38):     %d\n" % len(clin))
    both = set(curated) & set(clin)
    sys.stderr.write("in both:                     %d\n" % len(both))
    sys.stderr.write("ClinVar-only (new):          %d\n" % len(set(clin) - set(curated)))
    sys.stderr.write("curated-only (not in ClinVar): %d\n" % len(set(curated) - set(clin)))

    # Column names match variant_lookup.csv so variants.csv is a drop-in replacement
    # for every consumer (annotate_variants.py, vep_annotate.py) with no code changes.
    fields = ["HVGS", "Common name", "Gene", "Cohort_n", "Cohort frequency (%)",
              "Functionality", "ClinVar classification", "In Ithagenes? (Yes/No)",
              "In gnomAD? (Yes/No)", "gnomAD total AF (%)", "gnomAD ME AF (%)",
              "clinvar_review_status", "clinvar_id", "chrom", "pos", "ref", "alt", "source"]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for key in sorted(set(curated) | set(clin)):
            cu = curated.get(key, {})
            cv = clin.get(key, (None,))[0]
            srcs = [s for s, present in (("curated", cu), ("ClinVar", cv)) if present]
            w.writerow({
                "HVGS": key,
                "Common name": cu.get("Common name", ""),
                "Gene": cu.get("Gene") or (clin[key][2] if key in clin else ""),
                "Cohort_n": cu.get("Cohort_n", ""),
                "Cohort frequency (%)": cu.get("Cohort frequency (%)", ""),
                "Functionality": cu.get("Functionality", ""),
                "ClinVar classification": (cv.get("ClinicalSignificance", "") if cv else cu.get("ClinVar classification", "")),
                "In Ithagenes? (Yes/No)": cu.get("In Ithagenes? (Yes/No)", ""),
                "In gnomAD? (Yes/No)": cu.get("In gnomAD? (Yes/No)", ""),
                "gnomAD total AF (%)": cu.get("gnomAD total AF (%)", ""),
                "gnomAD ME AF (%)": cu.get("gnomAD ME AF (%)", ""),
                "clinvar_review_status": (cv.get("ReviewStatus", "") if cv else ""),
                "clinvar_id": (cv.get("VariationID", "") if cv else ""),
                "chrom": (cv.get("Chromosome", "") if cv else ""),
                "pos": (cv.get("PositionVCF", "") if cv else ""),
                "ref": (cv.get("ReferenceAlleleVCF", "") if cv else ""),
                "alt": (cv.get("AlternateAlleleVCF", "") if cv else ""),
                "source": "+".join(srcs),
            })
    sys.stderr.write("wrote %s\n" % a.out)


if __name__ == "__main__":
    main()
    