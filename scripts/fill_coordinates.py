#!/usr/bin/env python3
"""
fill_coordinates.py -- fill missing GRCh38 coordinates in naming_layer.csv using
VariantValidator, split into the two validated paths:

  * SUBSTITUTIONS (725): queried via the REST API using the MANE transcript
    c. notation. Coordinates extracted either from the direct grch38 VCF block
    or, for out-of-transcript promoter variants, parsed from the genomic
    re-submit suggestion in validation_warnings. This path was validated against
    123 known coordinates (45 general + 78 HBB promoter) with ZERO mismatches.

  * INDELS (361): the REST API's transcript endpoint returns only an HGVS-g range
    (e.g. g.5227099_5227100del) for out-of-transcript indels, whose coordinate
    convention differs from VCF (anchor-base / left-alignment). Rather than
    convert conventions by hand (error-prone), these are exported to a batch file
    for VariantValidator's WEB batch tool, which returns proper VCF coordinates.
    Import the returned coordinates with --merge-indels.

Usage:
  # step 1: fill substitutions + write the indel batch file
  python fill_coordinates.py --fill-subs

  # step 2 (after running the indel batch on variantvalidator.org):
  python fill_coordinates.py --merge-indels vv_batch_results.txt
"""
import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request

NAMING = "databases/naming_layer.csv"
GENE_TX = {"HBB": "NM_000518.5", "HBA1": "NM_000558.5", "HBA2": "NM_000517.6"}
API = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{v}/mane"


def is_sub(hgvs):
    return (">" in hgvs) and ("_" not in hgvs)


def is_out_of_transcript(hgvs):
    """Positions outside the coding transcript need the warning/two-step path:
    negative c. (5'UTR/promoter) AND c.* (3'UTR) positions."""
    return bool(re.search(r"c\.\[?\s*-\d", hgvs)) or bool(re.search(r"c\.\*\d", hgvs))


def is_compound(hgvs):
    return "[" in hgvs and ";" in hgvs


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
    """Direct VCF block (checked across ALL keys first), else substitution-warning
    genomic coordinate. Two passes so an empty response key (which occurs e.g.
    with LRG-pending warnings) cannot short-circuit a populated one."""
    # pass 1: any key with a direct GRCh38 VCF block
    for key, val in vv.items():
        if key.startswith("__") or not isinstance(val, dict):
            continue
        vcf = val.get("primary_assembly_loci", {}).get("grch38", {}).get("vcf", {})
        if vcf and vcf.get("pos"):
            return (str(vcf.get("chr", "")).replace("chr", ""),
                    str(vcf.get("pos", "")), vcf.get("ref", ""), vcf.get("alt", ""))
    # pass 2: substitution re-submit suggestion in any warning
    for key, val in vv.items():
        if key.startswith("__") or not isinstance(val, dict):
            continue
        for w in val.get("validation_warnings", []):
            m = re.search(r"NC_0000(\d+)\.\d+:g\.(\d+)([ACGT]+)>([ACGT]+)", w)
            if m:
                return (str(int(m.group(1))), m.group(2), m.group(3), m.group(4))
    return None


def hgvs_to_tx(hgvs):
    gene, rest = hgvs.split(":", 1)
    tx = GENE_TX.get(gene.strip())
    return f"{tx}:{rest.strip()}" if tx else None


def normalise_hgvs(hgvs):
    """Fix trivial, unambiguous HGVS typos before querying.
    Currently: a missing dot after 'c' (e.g. 'c224A>G' -> 'c.224A>G')."""
    # gene:cNNN... -> gene:c.NNN...
    return re.sub(r":c(?=[\d*\-])", ":c.", hgvs)


def clean_hgvs(hgvs):
    """Normalise indel HGVS quirks (validated: reproduces 30/30 known indel
    coordinates). delins is left intact (its inserted sequence is required)."""
    hgvs = normalise_hgvs(hgvs)
    if "delins" in hgvs:
        return hgvs
    # drop redundant deleted sequence and size suffixes: delACGT / delNNbp / delNN -> del
    hgvs = re.sub(r"(del)[ACGT]+$", r"\1", hgvs)
    hgvs = re.sub(r"(del)\d+bp$", r"\1", hgvs)
    hgvs = re.sub(r"(del)\d+$", r"\1", hgvs)
    return hgvs


def extract_genomic_hgvs_from_warning(vv):
    """For out-of-transcript variants, VV suggests a genomic HGVS in a warning:
    'Instead re-submit NC_0000XX.XX:g....'. Return that NC_ HGVS string."""
    for key, val in vv.items():
        if key.startswith("__") or not isinstance(val, dict):
            continue
        for w in val.get("validation_warnings", []):
            m = re.search(r"(NC_0000\d+\.\d+:g\.[\w>_]+)", w)
            if m:
                return m.group(1)
    return None


def query_genomic_vcf(nc_hgvs):
    """Query VV with a genomic NC_ HGVS and return the GRCh38 VCF tuple."""
    vv = query_vv(nc_hgvs)  # VV accepts genomic HGVS directly
    got = extract(vv)
    return got


def resolve_one(hgvs):
    """Resolve a single (non-compound) HGVS to (chrom,pos,ref,alt) via the
    validated paths: direct/warning for subs, direct for in-transcript indels,
    two-step for out-of-transcript indels. Returns None if unresolved."""
    tx = hgvs_to_tx(clean_hgvs(hgvs))
    if not tx:
        return None
    vv = query_vv(tx)
    got = extract(vv)
    if got and got[0] and got[1]:
        return got
    nc = extract_genomic_hgvs_from_warning(vv)
    if nc:
        got = query_genomic_vcf(nc)
        if got and got[0] and got[1]:
            return got
    return None


def fill_compound_alleles(pause):
    """Split gene:c.[a;b;c] into components, resolve each, store joined coords.
    Consumers key on HVGS, so joined coordinates (a;b;c) do not break them and
    honestly represent the compound nature. Protein-notation compounds (p.[...])
    are left for flagging (ambiguous at nucleotide level)."""
    rows = list(csv.DictReader(open(NAMING)))
    fields = list(rows[0].keys())
    targets = [r for r in rows if not r.get("chrom", "").strip()
               and is_compound(r["HGVS"]) and ":p." not in r["HGVS"]]
    print(f"compound alleles to split & fill: {len(targets)}\n")

    filled = partial = failed = 0
    fail_list = []
    for r in targets:
        hgvs = r["HGVS"].strip()
        gene = hgvs.split(":")[0]
        # extract the [a;b;c] payload
        m = re.search(r"c\.\[(.+)\]", hgvs)
        if not m:
            failed += 1; fail_list.append(hgvs); continue
        parts = [p.strip() for p in m.group(1).split(";")]
        coords = []
        for p in parts:
            comp = f"{gene}:c.{p}"
            got = resolve_one(comp)
            time.sleep(pause)
            coords.append(got)
        good = [c for c in coords if c]
        if len(good) == len(parts):
            r["chrom"] = ";".join(c[0] for c in coords)
            r["pos"] = ";".join(c[1] for c in coords)
            r["ref"] = ";".join(c[2] for c in coords)
            r["alt"] = ";".join(c[3] for c in coords)
            filled += 1
        elif good:
            # store the components we could resolve, mark partial
            r["chrom"] = ";".join(c[0] if c else "?" for c in coords)
            r["pos"] = ";".join(c[1] if c else "?" for c in coords)
            r["ref"] = ";".join(c[2] if c else "?" for c in coords)
            r["alt"] = ";".join(c[3] if c else "?" for c in coords)
            partial += 1
        else:
            failed += 1; fail_list.append(hgvs)

    with open(NAMING, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print(f"Compound alleles: {filled} fully filled, {partial} partial, "
          f"{failed} failed")
    for h in fail_list:
        print(f"    {h}")


def fill_out_of_transcript_indels(pause):
    """Two-step fill for out-of-transcript indels (negative c.-):
       transcript query -> genomic HGVS from warning -> genomic query -> VCF.
       Validated: reproduces known coordinate for HBB:c.-77_-76del (5227096)."""
    rows = list(csv.DictReader(open(NAMING)))
    fields = list(rows[0].keys())
    targets = [r for r in rows if not r.get("chrom", "").strip()
               and not is_sub(r.get("HGVS", ""))
               and is_out_of_transcript(r["HGVS"])
               and not is_compound(r["HGVS"])]
    print(f"out-of-transcript indels to fill: {len(targets)}\n")

    filled = failed = 0
    fail_list = []
    for r in targets:
        hgvs = clean_hgvs(r["HGVS"].strip())
        tx_hgvs = hgvs_to_tx(hgvs)
        if not tx_hgvs:
            failed += 1; fail_list.append(r["HGVS"]); continue
        vv = query_vv(tx_hgvs)
        time.sleep(pause)
        # try direct VCF first (some out-of-transcript indels return it directly)
        got = extract(vv)
        if not (got and got[0] and got[1]):
            # else two-step: genomic HGVS from warning -> genomic query -> VCF
            nc = extract_genomic_hgvs_from_warning(vv)
            if not nc:
                failed += 1; fail_list.append(r["HGVS"]); continue
            got = query_genomic_vcf(nc)
            time.sleep(pause)
        if got and got[0] and got[1]:
            r["chrom"], r["pos"], r["ref"], r["alt"] = got
            filled += 1
        else:
            failed += 1; fail_list.append(r["HGVS"])

    with open(NAMING, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print(f"Out-of-transcript indels: {filled} filled, {failed} failed")
    for h in fail_list:
        print(f"    {h}")


def fill_indels(pause):
    """Fill in-transcript indels via the API direct VCF block. Reports
    out-of-transcript (negative) and compound indels separately for other paths."""
    rows = list(csv.DictReader(open(NAMING)))
    fields = list(rows[0].keys())
    missing_indels = [r for r in rows
                      if not r.get("chrom", "").strip()
                      and not is_sub(r.get("HGVS", ""))]

    in_tx = [r for r in missing_indels
             if not is_out_of_transcript(r["HGVS"]) and not is_compound(r["HGVS"])]
    out_tx = [r for r in missing_indels if is_out_of_transcript(r["HGVS"])]
    compound = [r for r in missing_indels if is_compound(r["HGVS"])]

    print(f"indels missing coords: {len(missing_indels)}")
    print(f"  in-transcript (API): {len(in_tx)}  "
          f"out-of-transcript: {len(out_tx)}  compound: {len(compound)}\n")

    filled = failed = 0
    fail_list = []
    for i, r in enumerate(in_tx, 1):
        hgvs = clean_hgvs(r["HGVS"].strip())
        tx_hgvs = hgvs_to_tx(hgvs)
        if not tx_hgvs:
            failed += 1
            fail_list.append(r["HGVS"])
            continue
        got = extract(query_vv(tx_hgvs))
        time.sleep(pause)
        if got and got[0] and got[1]:
            r["chrom"], r["pos"], r["ref"], r["alt"] = got
            filled += 1
        else:
            failed += 1
            fail_list.append(r["HGVS"])
        if i % 50 == 0:
            print(f"  ...{i}/{len(in_tx)} ({filled} filled, {failed} failed)")

    with open(NAMING, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"\nIn-transcript indels: {filled} filled, {failed} failed")
    if fail_list:
        print("  failed (review):")
        for h in fail_list[:20]:
            print(f"    {h}")
    print(f"\nRemaining for other paths: {len(out_tx)} out-of-transcript, "
          f"{len(compound)} compound (handled separately)")


def fill_subs(pause):
    rows = list(csv.DictReader(open(NAMING)))
    fields = rows[0].keys()
    missing_subs = [r for r in rows
                    if not r.get("chrom", "").strip() and is_sub(r.get("HGVS", ""))]
    print(f"Filling {len(missing_subs)} substitution coordinates via VariantValidator...")

    filled = failed = 0
    fail_list = []
    for i, r in enumerate(missing_subs, 1):
        hgvs = normalise_hgvs(r["HGVS"].strip())
        tx_hgvs = hgvs_to_tx(hgvs)
        if not tx_hgvs:
            failed += 1
            fail_list.append(hgvs)
            continue
        got = extract(query_vv(tx_hgvs))
        time.sleep(pause)
        if got and got[0] and got[1]:
            r["chrom"], r["pos"], r["ref"], r["alt"] = got
            filled += 1
        else:
            failed += 1
            fail_list.append(hgvs)
        if i % 50 == 0:
            print(f"  ...{i}/{len(missing_subs)} ({filled} filled, {failed} failed)")

    # write back
    with open(NAMING, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fields))
        w.writeheader()
        w.writerows(rows)
    print(f"\nSubstitutions: {filled} filled, {failed} failed")
    if fail_list:
        print("  failed (left blank, review):")
        for h in fail_list[:20]:
            print(f"    {h}")


def write_indel_batch(path):
    rows = list(csv.DictReader(open(NAMING)))
    indels = [r for r in rows
              if not r.get("chrom", "").strip() and not is_sub(r.get("HGVS", ""))]
    with open(path, "w") as f:
        for r in indels:
            tx = hgvs_to_tx(r["HGVS"].strip())
            if tx:
                f.write(tx + "\n")
    print(f"Wrote {len(indels)} indel HGVS strings to {path}")
    print("  -> paste this file's contents into the VariantValidator web batch tool")
    print("     (https://variantvalidator.org, Batch tool), select GRCh38,")
    print("     download results, then run: --merge-indels <results>")


def merge_indels(results_path, pause):
    """Merge VCF coordinates from a VariantValidator batch-results file.
    Expects lines containing the transcript HGVS and a GRCh38 VCF coordinate;
    this parser is tolerant and matches by the NM_...:c. string."""
    rows = list(csv.DictReader(open(NAMING)))
    fields = list(rows[0].keys())
    # index naming rows by their transcript-form HGVS
    tx_index = {}
    for r in rows:
        if not r.get("chrom", "").strip():
            tx = hgvs_to_tx(r["HGVS"].strip())
            if tx:
                tx_index.setdefault(tx, []).append(r)

    filled = 0
    with open(results_path) as f:
        text = f.read()
    # find "NM_...:c.<...>" ... "chr-pos-ref-alt" or VCF fields nearby
    # (tolerant: match a transcript HGVS then the first GRCh38 VCF-like tuple after it)
    for line in text.splitlines():
        mtx = re.search(r"(NM_\d+\.\d+:c\.[^\s,\"]+)", line)
        mvcf = re.search(r"\b(\d{1,2})[-:](\d+)[-:]([ACGT]+)[-:]([ACGT]+)", line)
        if mtx and mvcf and mtx.group(1) in tx_index:
            for r in tx_index[mtx.group(1)]:
                r["chrom"], r["pos"], r["ref"], r["alt"] = (
                    mvcf.group(1), mvcf.group(2), mvcf.group(3), mvcf.group(4))
                filled += 1
    with open(NAMING, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"Merged {filled} indel coordinates from {results_path}")
    still = sum(1 for r in rows if not r.get("chrom", "").strip())
    print(f"Still missing coordinates: {still}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fill-subs", action="store_true")
    ap.add_argument("--fill-indels", action="store_true")
    ap.add_argument("--fill-outtx-indels", action="store_true")
    ap.add_argument("--fill-compound", action="store_true")
    ap.add_argument("--indel-batch", default="databases/indel_batch_for_vv.txt")
    ap.add_argument("--merge-indels", default="")
    ap.add_argument("--pause", type=float, default=0.4)
    args = ap.parse_args()

    if args.fill_subs:
        fill_subs(args.pause)
        write_indel_batch(args.indel_batch)
    elif args.fill_indels:
        fill_indels(args.pause)
    elif args.fill_outtx_indels:
        fill_out_of_transcript_indels(args.pause)
    elif args.fill_compound:
        fill_compound_alleles(args.pause)
    elif args.merge_indels:
        merge_indels(args.merge_indels, args.pause)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
