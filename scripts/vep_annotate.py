import requests
import csv
import sys
import time

input_vcf = sys.argv[1]
output_csv = sys.argv[2]
server = "https://rest.ensembl.org"
sample = input_vcf.split("/")[-1].split(".")[0]

import os
_cat_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                         "databases", "variants.csv")
CATALOGUE_HVGS = set()
try:
    with open(_cat_path, encoding="utf-8-sig") as _f:
        import csv as _csv
        for _row in _csv.DictReader(_f):
            _h = (_row.get("HVGS") or "").strip()
            if _h:
                CATALOGUE_HVGS.add(_h)
except Exception:
    pass

def normalise(hgvsc, gene):
    """ENST00000485743.2:c.-138C>A -> HBB:c.-138C>A"""
    if not hgvsc or ":c." not in hgvsc:
        return None
    c = hgvsc.split(":", 1)[1]
    return "%s:%s" % (gene, c)


def pick_consequence(transcript_consequences, catalogue_hvgs):
    globin = [tc for tc in transcript_consequences
              if tc.get("gene_symbol") in ("HBB", "HBA1", "HBA2")]
    if not globin:
        return None
    # 1. transcript whose GENE:c.notation is a known catalogue variant
    for tc in globin:
        key = normalise(tc.get("hgvsc", ""), tc["gene_symbol"])
        if key and key in catalogue_hvgs:
            return (tc["hgvsc"], tc["gene_symbol"],
                    tc.get("consequence_terms", ["unknown"])[0])
    # 2. canonical transcript with hgvsc
    for tc in globin:
        if tc.get("canonical") == 1 and tc.get("hgvsc"):
            return (tc["hgvsc"], tc["gene_symbol"],
                    tc.get("consequence_terms", ["unknown"])[0])
    # 3. any hgvsc
    for tc in globin:
        if tc.get("hgvsc"):
            return (tc["hgvsc"], tc["gene_symbol"],
                    tc.get("consequence_terms", ["unknown"])[0])
    # 4. canonical consequence, no hgvsc
    canon = next((tc for tc in globin if tc.get("canonical") == 1), globin[0])
    return ("unknown", canon["gene_symbol"],
            canon.get("consequence_terms", ["unknown"])[0])
            
def query_vep(chrom, pos, ref, alt):
    chrom_num = chrom.replace("chr", "")
    ref_len = len(ref)
    alt_len = len(alt)
    
    if ref_len == 1 and alt_len == 1:
        # SNV
        url = f"{server}/vep/human/region/{chrom_num}:{pos}:{pos}/{alt}"
    elif ref_len > alt_len:
        # Deletion
        del_start = int(pos) + len(alt)
        del_end = int(pos) + len(ref) - 1
        url = f"{server}/vep/human/region/{chrom_num}:{del_start}-{del_end}:1/-"
    elif alt_len > ref_len:
        # Insertion
        ins_seq = alt[len(ref):]
        ins_pos = int(pos) + len(ref) - 1
        url = f"{server}/vep/human/region/{chrom_num}:{ins_pos}-{ins_pos + 1}:1/{ins_seq}"
    else:
        # Complex, try as-is
        url = f"{server}/vep/human/region/{chrom_num}:{pos}:{pos}/{alt}"
    
    try:
        r = requests.get(url, headers={"Content-Type": "application/json"}, params={"hgvs": 1, "regulatory": 1, "canonical": 1})
        time.sleep(0.1)  # rate limit
        if r.ok:
            data = r.json()
            if data and "transcript_consequences" in data[0]:
                result = pick_consequence(data[0]["transcript_consequences"], CATALOGUE_HVGS)
                if result:
                    return result
            if data and "regulatory_feature_consequences" in data[0]:
                rfc = data[0]["regulatory_feature_consequences"][0]
                return (f"regulatory:{rfc.get('biotype','regulatory')}", "",
                        rfc.get("consequence_terms", ["regulatory_region_variant"])[0])
            return "non-globin", "", ""
        return "api_error", "", ""
    except Exception as e:
        return str(e), "", ""

with open(input_vcf) as f, open(output_csv, 'w') as out:
    writer = csv.writer(out)
    writer.writerow(["Sample", "Chromosome", "Position", "Ref", "Alt", "HGVS", "Gene", "Consequence", "Genotype", "Quality"])
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        chrom, pos, ref, alt, qual = parts[0], parts[1], parts[3], parts[4], parts[5]
        gt = parts[9].split(":")[0]
        hgvs, gene, consequence = query_vep(chrom, pos, ref, alt)
        writer.writerow([sample, chrom, pos, ref, alt, hgvs, gene, consequence, gt, qual])
        

            