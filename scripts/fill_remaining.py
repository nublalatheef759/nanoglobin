#!/usr/bin/env python3
"""
fill_remaining.py -- resolve the last fillable flagged variants:

  * boundary-spanning deletions whose only problem was a BACKWARDS position range
    (VariantValidator reports 'end < start' and suggests the reordered form).
    We query VV, read its LovdSyntaxcheckSuggestions reordered HGVS, and re-query
    that -- letting VariantValidator itself supply the corrected description.

Leaves genuinely-unresolvable flags (protein notation, dual-gene uncertainty)
untouched.
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request

NAMING = "databases/naming_layer.csv"
GENE_TX = {"HBB": "NM_000518.5", "HBA1": "NM_000558.5", "HBA2": "NM_000517.6"}
API = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{v}/mane"


def query_vv(variant):
    url = API.format(v=urllib.parse.quote(variant, safe=""))
    try:
        req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.load(resp)
    except Exception as e:
        return {"__error__": str(e)}


def extract(vv):
    for key, val in vv.items():
        if key.startswith("__") or not isinstance(val, dict):
            continue
        vcf = val.get("primary_assembly_loci", {}).get("grch38", {}).get("vcf", {})
        if vcf and vcf.get("pos"):
            return (str(vcf.get("chr", "")).replace("chr", ""),
                    str(vcf.get("pos", "")), vcf.get("ref", ""), vcf.get("alt", ""))
    return None


def get_suggestion(vv):
    """Read VariantValidator's LovdSyntaxcheckSuggestions reordered HGVS."""
    for key, val in vv.items():
        if key.startswith("__") or not isinstance(val, dict):
            continue
        for w in val.get("validation_warnings", []):
            m = re.search(r"suggestion = (NM_\d+\.\d+:c\.[^\s,\]]+)", w)
            if m:
                return m.group(1)
    return None


def main():
    rows = list(csv.DictReader(open(NAMING)))
    fields = list(rows[0].keys())

    # target: rows still missing coords whose flag is a boundary-spanning deletion
    targets = [r for r in rows if not r.get("chrom", "").strip()
               and "boundary_spanning" in r.get("coord_status", "")]
    print(f"boundary-spanning deletions to retry: {len(targets)}\n")

    filled = 0
    for r in targets:
        hgvs = r["HGVS"].strip()
        gene = hgvs.split(":")[0]
        if gene not in GENE_TX:
            continue
        tx = f"{GENE_TX[gene]}:{hgvs.split(':', 1)[1]}"
        vv = query_vv(tx)
        time.sleep(0.4)
        got = extract(vv)
        if not got:
            # try VV's own reordered suggestion
            sug = get_suggestion(vv)
            if sug:
                got = extract(query_vv(sug))
                time.sleep(0.4)
        if got and got[0] and got[1]:
            r["chrom"], r["pos"], r["ref"], r["alt"] = got
            r["coord_status"] = ""  # resolved
            filled += 1
            print(f"  FILLED {hgvs} -> {got[0]}:{got[1]}")
        else:
            print(f"  still unresolved: {hgvs}")

    with open(NAMING, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    wc = sum(1 for r in rows if r.get("chrom", "").strip())
    print(f"\nfilled {filled} boundary-spanning deletions")
    print(f"coverage: {wc}/{len(rows)} ({100*wc*10//len(rows)/10}%), "
          f"missing: {len(rows) - wc}")


if __name__ == "__main__":
    main()
