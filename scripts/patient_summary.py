import csv
import sys
import re

input_file = sys.argv[1]
output_file = sys.argv[2]

patients = {}

with open(input_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = row["Sample"]
        if sample not in patients:
            patients[sample] = {"HBA": set(), "HBB": set(), "HBA_zyg": set(), "HBB_zyg": set()}
        
        if row["Tool"] == "Coverage":
            continue
        if row["Consequence"] == "upstream_gene_variant":
            continue
        
        gt = row["Genotype"]
        if gt == "1/1":
            zyg = "homozygous"
        elif gt == "0/1":
            zyg = "heterozygous"
        else:
            zyg = "uncertain"
        
        hgvs = row["HGVS"]
        if not hgvs or hgvs in ["unknown", "non-globin", ""]:
            if "SV_" not in row.get("Consequence", ""):
                continue
            hgvs = row["HGVS"] if row["HGVS"] else row["Consequence"]
        
        # Strip overlap percentage for dedup
        base_name = re.sub(r'\s*\(\d+% overlap\)', '', hgvs)
        # Skip unknowns
        if base_name.startswith("unknown"):
            continue
        
        gene = row["Gene"]
        if gene in ["HBA", "HBA1", "HBA2"]:
            patients[sample]["HBA"].add(base_name)
            patients[sample]["HBA_zyg"].add(zyg)
        elif gene == "HBB":
            patients[sample]["HBB"].add(base_name)
            patients[sample]["HBB_zyg"].add(zyg)

with open(output_file, 'w') as out:
    writer = csv.writer(out)
    writer.writerow(["Sample", "HBA", "HBA_zygosity", "HBB", "HBB_zygosity"])
    for sample in sorted(patients):
        p = patients[sample]
        writer.writerow([
            sample,
            "; ".join(sorted(p["HBA"])) if p["HBA"] else "wild-type",
            "; ".join(sorted(p["HBA_zyg"])) if p["HBA_zyg"] else "",
            "; ".join(sorted(p["HBB"])) if p["HBB"] else "wild-type",
            "; ".join(sorted(p["HBB_zyg"])) if p["HBB_zyg"] else ""
        ])
