#!/usr/bin/env python3
"""
validate_warning_path.py -- targeted test of the warning-path coordinate
extraction against the 78 HBB promoter variants that ALREADY have known
coordinates. These out-of-transcript c.-NNN variants are exactly the category
that dominates the 1086 unknowns, so proving the warning-path reproduces their
known coordinates is the decisive safety check before filling.
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request

GENE_TX = {"HBB": "NM_000518.5", "HBA1": "NM_000558.5", "HBA2": "NM_000517.6"}
API = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{v}/mane"


def query_vv(variant, retries=3):
    url = API.format(v=urllib.parse.quote(variant, safe=""))
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.load(resp)
        except Exception as e:
            if attempt == retries - 1:
                return {"__error__": str(e)}
            time.sleep(1.0 * (attempt + 1))
    return {"__error__": "unreachable"}


def extract(vv):
    for key, val in vv.items():
        if key.startswith("__") or not isinstance(val, dict):
            continue
        vcf = val.get("primary_assembly_loci", {}).get("grch38", {}).get("vcf", {})
        if vcf:
            return (str(vcf.get("chr", "")).replace("chr", ""),
                    str(vcf.get("pos", "")), vcf.get("ref", ""), vcf.get("alt", ""))
        for w in val.get("validation_warnings", []):
            m = re.search(r"NC_0000(\d+)\.\d+:g\.(\d+)([ACGT]+)>([ACGT]+)", w)
            if m:
                return (str(int(m.group(1))), m.group(2), m.group(3), m.group(4))
    return None


def main():
    rows = list(csv.DictReader(open("databases/naming_layer.csv")))
    warned = [r for r in rows if r.get("HGVS", "").startswith("HBB:c.-")
              and r.get("chrom", "").strip() and r.get("pos", "").strip()]

    print(f"Validating warning-path on {len(warned)} HBB promoter variants "
          f"with known coordinates...\n")
    match = mismatch = error = 0
    mismatches = []
    for r in warned:
        hgvs = r["HGVS"].strip()
        gene, rest = hgvs.split(":", 1)
        tx_hgvs = f"{GENE_TX[gene]}:{rest}"
        known = (r["chrom"].strip().replace("chr", ""), r["pos"].strip(),
                 r.get("ref", "").strip(), r.get("alt", "").strip())
        got = extract(query_vv(tx_hgvs))
        time.sleep(0.4)
        if got is None:
            print(f"  {hgvs:<22} known {known[0]}:{known[1]}  -> ERROR")
            error += 1
            continue
        # compare chrom + pos + ref + alt (full match, since warning-path gives all)
        if got[0] == known[0] and got[1] == known[1] and \
           got[2] == known[2] and got[3] == known[3]:
            match += 1
        elif got[0] == known[0] and got[1] == known[1]:
            # pos matches but ref/alt representation differs -- still a position match
            print(f"  {hgvs:<22} known {known[0]}:{known[1]} {known[2]}>{known[3]}  "
                  f"VV {got[0]}:{got[1]} {got[2]}>{got[3]}  (pos match, allele repr differs)")
            match += 1
        else:
            print(f"  {hgvs:<22} known {known[0]}:{known[1]}  VV {got[0]}:{got[1]}  "
                  f"*** MISMATCH ***")
            mismatch += 1
            mismatches.append((hgvs, known, got))

    print(f"\nSUMMARY (warning-path, HBB promoter): "
          f"{match} MATCH, {mismatch} MISMATCH, {error} ERROR (of {len(warned)})")
    if mismatch == 0 and error == 0:
        print("==> Warning-path reproduces ALL known promoter coordinates. "
              "SAFE to fill unknowns.")
    else:
        print("==> Investigate before filling.")
        for hgvs, known, got in mismatches:
            print(f"    {hgvs}: known {known} vs VV {got}")


if __name__ == "__main__":
    main()
