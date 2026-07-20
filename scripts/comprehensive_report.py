import csv
import glob
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.identify_sv import identify_sv

output_file = sys.argv[1]

with open(output_file, 'w') as out:
    writer = csv.writer(out)
    writer.writerow(["Sample", "Tool", "Chromosome", "Position", "Ref", "Alt", "HGVS", "Gene", "Consequence", "Genotype", "Quality"])

    for f in sorted(glob.glob("results/*.annotated.csv")):
        with open(f) as ann:
            reader = csv.reader(ann)
            next(reader)
            for row in reader:
                if row[5] == "non-globin":
                    continue
                writer.writerow([row[0], "Clair3", row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]])

    for f in sorted(glob.glob("variants/*/*.sniffles.vcf")):
        sample = f.split("/")[1]
        with open(f) as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                chrom, pos = parts[0], int(parts[1])
                in_hbb = chrom == "chr11" and 5225000 <= pos <= 5228000
                in_hba = chrom == "chr16" and 165000 <= pos <= 180000
                if not (in_hbb or in_hba):
                    continue
                svtype = ""
                svlen = ""
                for info in parts[7].split(";"):
                    if info.startswith("SVTYPE="):
                        svtype = info.split("=")[1]
                    if info.startswith("SVLEN="):
                        svlen = info.split("=")[1]
                if svtype in ["INV", "BND"]:
                    continue
                gt = parts[9].split(":")[0] if len(parts) > 9 else ""
                gene = "HBB" if in_hbb else "HBA"
                sv_name = identify_sv(chrom, str(pos), svtype, svlen)
                ref_display = parts[3][:20] + "..." if len(parts[3]) > 20 else parts[3]
                writer.writerow([sample, "Sniffles", chrom, pos, ref_display, svtype, sv_name, gene, "SV_" + svtype + "_" + svlen + "bp", gt, parts[5]])

    for f in sorted(glob.glob("variants/*/*.cutesv.vcf")):
        sample = f.split("/")[1]
        with open(f) as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                chrom, pos = parts[0], int(parts[1])
                in_hbb = chrom == "chr11" and 5225000 <= pos <= 5228000
                in_hba = chrom == "chr16" and 170000 <= pos <= 178000
                if not (in_hbb or in_hba):
                    continue
                svtype = ""
                svlen = ""
                for info in parts[7].split(";"):
                    if info.startswith("SVTYPE="):
                        svtype = info.split("=")[1]
                    if info.startswith("SVLEN="):
                        svlen = info.split("=")[1]
                if svtype in ["INV", "BND"]:
                    continue
                gt = parts[9].split(":")[0] if len(parts) > 9 else ""
                gene = "HBB" if in_hbb else "HBA"
                sv_name = identify_sv(chrom, str(pos), svtype, svlen)
                ref_display = parts[3][:20] + "..." if len(parts[3]) > 20 else parts[3]
                writer.writerow([sample, "CuteSV", chrom, pos, ref_display, svtype, sv_name, gene, "SV_" + svtype + "_" + svlen + "bp", gt, parts[5]])

    for f in sorted(glob.glob("variants/*/*.coverage.tsv")):
        sample = f.split("/")[1]
        with open(f) as cov:
            for line in cov:
                parts = line.strip().split("\t")
                writer.writerow([sample, "Coverage", "", "", "", "", "", parts[0], parts[1], "", ""])
