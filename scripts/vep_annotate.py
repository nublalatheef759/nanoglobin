import requests
import csv
import sys
import time

input_vcf = sys.argv[1]
output_csv = sys.argv[2]
server = "https://rest.ensembl.org"
sample = input_vcf.split("/")[-1].split(".")[0]

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
        r = requests.get(url, headers={"Content-Type": "application/json"}, params={"hgvs": 1})
        time.sleep(0.1)  # rate limit
        if r.ok:
            data = r.json()
            if data and "transcript_consequences" in data[0]:
                for tc in data[0]["transcript_consequences"]:
                    if tc.get("gene_symbol") in ["HBB", "HBA1", "HBA2"]:
                        return tc.get("hgvsc", "unknown"), tc.get("gene_symbol", ""), tc.get("consequence_terms", ["unknown"])[0]
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

            