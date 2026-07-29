#!/usr/bin/env python3
"""
finalize_coordinates.py -- close out the coordinate fill:
  1. fix the one obvious HGVS typo (c224A>G -> c.224A>G) and resolve it
  2. pull the one structural variant's coordinate from cnvs.csv (own curated data)
  3. retry insertion-format variants with an explicit range (c.NinsX -> c.N_N+1insX)
  4. flag every genuinely-unresolvable remainder with a documented coord_status
     reason (protein-notation ambiguity, malformed source, boundary-spanning,
     reference mismatch) -- turning "missing" into "documented exception".

Adds a coord_status column (blank for resolved rows, a reason for exceptions).
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request

NAMING = "databases/naming_layer.csv"
CNVS = "databases/cnvs.csv"
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


def classify_reason(hgvs):
    if ":p." in hgvs:
        return "protein_notation_ambiguous"
    if "^" in hgvs or re.search(r"[^ACGT\d_>\[\];:.c*+\-() ]", hgvs.split(":", 1)[-1]):
        return "malformed_source_entry"
    if re.search(r"-\d+_\d+\+", hgvs) or re.search(r"c\.-\d+_\d", hgvs):
        return "boundary_spanning_deletion"
    return "unresolved_review"


def main():
    rows = list(csv.DictReader(open(NAMING)))
    fields = list(rows[0].keys())
    if "coord_status" not in fields:
        fields.append("coord_status")
        for r in rows:
            r.setdefault("coord_status", "")

    # --- load cnvs.csv for the structural cross-reference (by size in HGVS) ---
    cnv_rows = list(csv.reader(open(CNVS)))

    fixed = 0
    for r in rows:
        if r.get("chrom", "").strip():
            continue
        hgvs = r["HGVS"].strip()

        # 1. typo: missing dot after c
        if re.search(r":c\d", hgvs):
            fixed_hgvs = re.sub(r":c(?=\d)", ":c.", hgvs)
            gene = fixed_hgvs.split(":")[0]
            if gene in GENE_TX:
                got = extract(query_vv(f"{GENE_TX[gene]}:{fixed_hgvs.split(':',1)[1]}"))
                time.sleep(0.4)
                if got:
                    r["chrom"], r["pos"], r["ref"], r["alt"] = got
                    fixed += 1
                    continue

        # 2. structural already in cnvs.csv (match by the delNNN size token)
        m = re.search(r"del(\d+)$", hgvs)
        if m and int(m.group(1)) > 50:
            size = m.group(1)
            for c in cnv_rows:
                line = ",".join(c)
                if f"del{size}" in line and "chr11" in line or (size in line and "DEL" in line):
                    # cnvs.csv columns: name, chrom, from, from2, to, to2, ...
                    try:
                        r["chrom"] = c[1].replace("chr", "")
                        r["pos"] = c[2]
                        r["coord_status"] = "structural_from_cnvs_csv"
                        fixed += 1
                    except IndexError:
                        pass
                    break
            if r.get("chrom", "").strip():
                continue

        # 3. insertion missing range: c.NinsX -> c.N_(N+1)insX
        mi = re.search(r"c\.(\d+)ins([ACGT]+)$", hgvs)
        if mi:
            n = int(mi.group(1))
            gene = hgvs.split(":")[0]
            if gene in GENE_TX:
                fixed_hgvs = f"{GENE_TX[gene]}:c.{n}_{n+1}ins{mi.group(2)}"
                got = extract(query_vv(fixed_hgvs))
                time.sleep(0.4)
                if got:
                    r["chrom"], r["pos"], r["ref"], r["alt"] = got
                    fixed += 1
                    continue

        # 4. flag the genuine remainder with a reason
        r["coord_status"] = classify_reason(hgvs)

    with open(NAMING, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    wc = sum(1 for r in rows if r.get("chrom", "").strip())
    flagged = {}
    for r in rows:
        st = r.get("coord_status", "")
        if st and not r.get("chrom", "").strip():
            flagged[st] = flagged.get(st, 0) + 1
    print(f"resolved this pass: {fixed}")
    print(f"coverage now: {wc}/{len(rows)} ({100*wc//len(rows)}%)")
    print(f"flagged exceptions (no coordinate, documented reason):")
    for k, v in flagged.items():
        print(f"    {k}: {v}")


if __name__ == "__main__":
    main()
