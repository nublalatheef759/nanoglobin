#!/usr/bin/env python3
"""
validate_indel_cleaning.py -- before using cleaned HGVS (dropping redundant
deleted bases, normalising insertions) to recover the ~8 indel failures, verify
the cleaning reproduces KNOWN indel coordinates. Tests a sample of the 385
known-coordinate indels: sends the cleaned HGVS to VariantValidator and checks
it returns the coordinate already on file. Same validate-before-trust gate used
for substitutions.
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request

GENE_TX = {"HBB": "NM_000518.5", "HBA1": "NM_000558.5", "HBA2": "NM_000517.6"}
API = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{v}/mane"


def clean_hgvs(hgvs):
    """Normalise indel HGVS quirks:
      - drop redundant deleted sequence: c.NNN_MMMdelACGT -> c.NNN_MMMdel
      - drop size suffixes: delNNbp / delNNN -> del
    (insertions and delins are left intact -- their sequence is required)
    """
    # delins must be preserved (has required inserted seq); only bare del is trimmed
    if "delins" in hgvs:
        return hgvs
    # c.<range>del<bases>  ->  c.<range>del   (bases redundant)
    hgvs = re.sub(r"(del)[ACGT]+$", r"\1", hgvs)
    # c.<range>del<NN>bp   ->  c.<range>del
    hgvs = re.sub(r"(del)\d+bp$", r"\1", hgvs)
    hgvs = re.sub(r"(del)\d+$", r"\1", hgvs)
    return hgvs


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


def main():
    rows = list(csv.DictReader(open("databases/naming_layer.csv")))
    known = [r for r in rows if r.get("chrom", "").strip()
             and ("del" in r.get("HGVS", "") or "ins" in r.get("HGVS", ""))
             and "_" in r.get("HGVS", "")]

    # prefer ones with redundant-del notation (the quirk we're validating)
    redundant = [r for r in known if re.search(r"del[ACGT]+$", r["HGVS"])]
    sample = (redundant[:20] + known[:10])[:30]

    print(f"Validating cleaned-HGVS indel path on {len(sample)} known indels...\n")
    match = mismatch = error = 0
    for r in sample:
        hgvs = r["HGVS"].strip()
        cleaned = clean_hgvs(hgvs)
        gene = hgvs.split(":")[0]
        if gene not in GENE_TX:
            continue
        tx = f"{GENE_TX[gene]}:{cleaned.split(':', 1)[1]}"
        known_c = (r["chrom"].strip(), r["pos"].strip())
        got = extract(query_vv(tx))
        time.sleep(0.4)
        if got is None:
            print(f"  {hgvs:<34} -> ERROR (cleaned: {cleaned})")
            error += 1
        elif got[0] == known_c[0] and got[1] == known_c[1]:
            match += 1
        else:
            print(f"  {hgvs:<34} known {known_c[0]}:{known_c[1]} "
                  f"vs VV {got[0]}:{got[1]}  *** MISMATCH ***")
            mismatch += 1

    print(f"\nSUMMARY: {match} MATCH, {mismatch} MISMATCH, {error} ERROR "
          f"(of {len(sample)})")
    if mismatch == 0 and error == 0:
        print("==> Cleaned-HGVS indel path reproduces known coordinates. SAFE.")
    else:
        print("==> Investigate before using cleaned-HGVS path.")


if __name__ == "__main__":
    main()
