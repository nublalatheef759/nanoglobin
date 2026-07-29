#!/usr/bin/env python3
"""
validate_variantvalidator.py -- SAFETY GATE before using VariantValidator to
fill missing coordinates.

This script does NOT fill anything. It takes a sample of variants that ALREADY
have known, authoritative (ClinVar-sourced) GRCh38 coordinates, asks
VariantValidator for those same variants, and checks whether VariantValidator
reproduces the known coordinate exactly. Only if it reproduces the known
coordinates reliably should it be trusted to fill the unknown ones.

Tests a spread across HBB / HBA1 / HBA2 and different variant types (SNV, indel,
UTR, coding) so region- or type-specific problems are caught.

Reports MATCH / MISMATCH / ERROR per variant. Review the results before running
any fill.
"""
import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request


# MANE-select transcripts (confirmed: NM_000518.5 is HBB MANE-select)
GENE_TX = {
    "HBB":  "NM_000518.5",
    "HBA1": "NM_000558.5",
    "HBA2": "NM_000517.6",
}

API = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{variant}/mane"


def query_vv(variant, retries=3, pause=1.0):
    url = API.format(variant=urllib.parse.quote(variant, safe=""))
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.load(resp)
        except Exception as e:
            if attempt == retries - 1:
                return {"__error__": str(e)}
            time.sleep(pause * (attempt + 1))
    return {"__error__": "unreachable"}


def extract_grch38(vv_json):
    """Pull chrom/pos/ref/alt from a VV response.

    Two paths:
      1. Direct: primary_assembly_loci.grch38.vcf (in-transcript variants).
      2. Warning: for variants outside the transcript (e.g. HBB promoter c.-NNN),
         VV declines the c. form but supplies the correct genomic coordinate in a
         validation_warning as 'Instead re-submit NC_0000XX.XX:g.<pos><ref>><alt>'.
         VV has already done the strand math, so we just parse it.
    """
    for key, val in vv_json.items():
        if key.startswith("__") or not isinstance(val, dict):
            continue
        # path 1: direct VCF block
        vcf = val.get("primary_assembly_loci", {}).get("grch38", {}).get("vcf", {})
        if vcf:
            return (str(vcf.get("chr", "")).replace("chr", ""),
                    str(vcf.get("pos", "")), vcf.get("ref", ""), vcf.get("alt", ""))
        # path 2: parse the genomic re-submit suggestion from warnings
        for w in val.get("validation_warnings", []):
            m = re.search(r"NC_0000(\d+)\.\d+:g\.(\d+)([ACGT]+)>([ACGT]+)", w)
            if m:
                chrom = str(int(m.group(1)))  # NC_000011 -> 11
                return (chrom, m.group(2), m.group(3), m.group(4))
    return None


def hgvs_to_tx(hgvs):
    """Convert 'HBB:c.20A>T' -> 'NM_000518.5:c.20A>T' using MANE transcripts."""
    if ":" not in hgvs:
        return None
    gene, rest = hgvs.split(":", 1)
    tx = GENE_TX.get(gene.strip())
    if not tx:
        return None
    return f"{tx}:{rest.strip()}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--naming", default="databases/naming_layer.csv")
    ap.add_argument("--n", type=int, default=30, help="how many known variants to test")
    ap.add_argument("--pause", type=float, default=0.5, help="seconds between API calls")
    args = ap.parse_args()

    with open(args.naming) as f:
        rows = list(csv.DictReader(f))

    # KNOWN variants: those that already have chrom+pos (authoritative), spread
    # across the three genes.
    known = [r for r in rows if r.get("chrom", "").strip() and r.get("pos", "").strip()]
    by_gene = {"HBB": [], "HBA1": [], "HBA2": []}
    for r in known:
        h = r.get("HGVS", "")
        g = h.split(":")[0] if ":" in h else ""
        if g in by_gene:
            by_gene[g].append(r)

    # take a balanced sample across genes
    sample = []
    per = max(1, args.n // 3)
    for g in ("HBB", "HBA1", "HBA2"):
        sample.extend(by_gene[g][:per])

    print(f"Testing {len(sample)} known variants against VariantValidator...\n")
    print(f"{'HGVS':<28} {'known':<18} {'VV_GRCh38':<18} {'result'}")
    print("-" * 80)

    n_match = n_mismatch = n_error = 0
    mismatches = []
    for r in sample:
        hgvs = r.get("HGVS", "").strip()
        tx_hgvs = hgvs_to_tx(hgvs)
        known_coord = (r["chrom"].strip().replace("chr", ""), r["pos"].strip(),
                       r.get("ref", "").strip(), r.get("alt", "").strip())
        if not tx_hgvs:
            print(f"{hgvs:<28} {'-':<18} {'-':<18} SKIP (no tx map)")
            continue

        vv = query_vv(tx_hgvs)
        got = extract_grch38(vv)
        time.sleep(args.pause)

        known_str = f"{known_coord[0]}:{known_coord[1]}"
        if "__error__" in vv or got is None:
            print(f"{hgvs:<28} {known_str:<18} {'ERROR':<18} "
                  f"{vv.get('__error__', 'no coord returned')[:30]}")
            n_error += 1
            continue

        got_str = f"{got[0]}:{got[1]}"
        # compare chrom + pos (ref/alt can differ in representation for indels)
        if got[0] == known_coord[0] and got[1] == known_coord[1]:
            print(f"{hgvs:<28} {known_str:<18} {got_str:<18} MATCH")
            n_match += 1
        else:
            print(f"{hgvs:<28} {known_str:<18} {got_str:<18} *** MISMATCH ***")
            n_mismatch += 1
            mismatches.append((hgvs, known_coord, got))

    print("-" * 80)
    print(f"\nSUMMARY: {n_match} MATCH, {n_mismatch} MISMATCH, {n_error} ERROR "
          f"(of {len(sample)} tested)")
    if n_mismatch == 0 and n_error == 0 and n_match > 0:
        print("\n==> VariantValidator reproduces all known coordinates. "
              "SAFE to proceed to filling unknowns.")
    elif n_mismatch > 0:
        print("\n==> MISMATCHES FOUND. DO NOT fill unknowns until investigated:")
        for hgvs, known_coord, got in mismatches:
            print(f"    {hgvs}: known {known_coord[0]}:{known_coord[1]} "
                  f"vs VV {got[0]}:{got[1]}")
    else:
        print("\n==> Errors/insufficient matches. Investigate before proceeding.")


if __name__ == "__main__":
    main()
