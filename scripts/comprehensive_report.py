import csv
import glob
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.identify_sv import identify_sv

output_file = sys.argv[1]
if len(sys.argv) < 3:
    sys.exit("usage: comprehensive_report.py <output.csv> <sample> [sample ...]")
SAMPLES = sys.argv[2:]

def _declared(pattern, extract):
    """Glob, then keep only files whose sample name is in the declared set."""
    keep = []
    for f in sorted(glob.glob(pattern)):
        s = extract(f)
        if s in SAMPLES:
            keep.append(f)
    missing = set(SAMPLES) - {extract(f) for f in keep}
    if missing:
        sys.exit("ERROR: no input found for declared sample(s): %s (%s)"
                 % (", ".join(sorted(missing)), pattern))
    return keep

_by_basename = lambda f: os.path.basename(f).split(".")[0]
_by_dirname  = lambda f: os.path.basename(os.path.dirname(f))

with open(output_file, 'w') as out:
    writer = csv.writer(out)
    writer.writerow(["Sample", "Tool", "Chromosome", "Position", "Ref", "Alt", "HGVS", "Gene", "Consequence", "Genotype", "Quality"])

    for f in _declared("results/*.annotated.csv", _by_basename):
        with open(f) as ann:
            reader = csv.reader(ann)
            next(reader)
            for row in reader:
                if row[5] == "non-globin":
                    continue
                writer.writerow([row[0], "Clair3", row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]])

    for f in _declared("variants/*/*.sniffles.vcf", _by_dirname):
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

    for f in _declared("variants/*/*.cutesv.vcf", _by_dirname):
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

    for f in _declared("variants/*/*.coverage.tsv", _by_dirname):
        sample = f.split("/")[1]
        with open(f) as cov:
            for line in cov:
                parts = line.strip().split("\t")
                writer.writerow([sample, "Coverage", "", "", "", "", "", parts[0], parts[1], "", ""])
