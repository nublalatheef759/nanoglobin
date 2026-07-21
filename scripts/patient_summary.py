import csv
import sys
import re

input_file = sys.argv[1]
output_file = sys.argv[2]

# Rank zygosity calls so that, when multiple callers report the SAME variant
# (e.g. Sniffles 0/1 and CuteSV ./. for one -a3.7 deletion), the confident call
# wins and the missing-genotype one does not produce a spurious 'undetermined'.
ZYG_RANK = {"homozygous": 3, "heterozygous": 3, "reference": 2, "undetermined": 1}


def zygosity_from_gt(gt):
    g = (gt or "").strip().replace("|", "/")
    alleles = g.split("/")
    if alleles == ["1", "1"]:
        return "homozygous"
    if alleles in (["0", "1"], ["1", "0"]):
        return "heterozygous"
    if alleles == ["0", "0"]:
        return "reference"
    return "undetermined"


# per sample -> per gene -> {variant_name: best_zygosity}
patients = {}

with open(input_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = row["Sample"]
        if sample not in patients:
            patients[sample] = {"HBA": {}, "HBB": {}}

        if row["Tool"] == "Coverage":
            continue
        if row["Consequence"] == "upstream_gene_variant":
            continue

        zyg = zygosity_from_gt(row["Genotype"])

        hgvs = row["HGVS"]
        if not hgvs or hgvs in ["unknown", "non-globin", ""]:
            if "SV_" not in row.get("Consequence", ""):
                continue
            hgvs = row["HGVS"] if row["HGVS"] else row["Consequence"]

        # Strip overlap percentage for dedup
        base_name = re.sub(r'\s*\(\d+% overlap\)', '', hgvs)
        if base_name.startswith("unknown"):
            continue

        gene = row["Gene"]
        if gene in ["HBA", "HBA1", "HBA2"]:
            key = "HBA"
        elif gene == "HBB":
            key = "HBB"
        else:
            continue

        # Same variant may be reported by >1 caller: keep the highest-ranked
        # zygosity so a confident 0/1 is not overwritten by a ./. undetermined.
        variants = patients[sample][key]
        if base_name not in variants or ZYG_RANK[zyg] > ZYG_RANK[variants[base_name]]:
            variants[base_name] = zyg


with open(output_file, 'w') as out:
    writer = csv.writer(out)
    writer.writerow(["Sample", "HBA", "HBA_zygosity", "HBB", "HBB_zygosity"])
    for sample in sorted(patients):
        hba = patients[sample]["HBA"]
        hbb = patients[sample]["HBB"]
        hba_names = "; ".join(sorted(hba)) if hba else "wild-type"
        hbb_names = "; ".join(sorted(hbb)) if hbb else "wild-type"
        # zygosities listed in the SAME order as the names they belong to
        hba_zyg = "; ".join(hba[n] for n in sorted(hba)) if hba else ""
        hbb_zyg = "; ".join(hbb[n] for n in sorted(hbb)) if hbb else ""
        writer.writerow([sample, hba_names, hba_zyg, hbb_names, hbb_zyg])
        
