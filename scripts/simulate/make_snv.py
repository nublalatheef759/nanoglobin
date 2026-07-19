#!/usr/bin/env python3
"""
make_snv.py -- introduce a single point mutation into a simulation reference and
emit a truth VCF. Companion to make_haplotype.py (which does deletions).

Coordinates are given in FULL reference space; the contig offset embedded in the
sim_target.fa header (e.g. 'chr11:5220000-5232000') is handled internally.

Usage:
  python make_snv.py --ref sim_target.fa --chrom chr11 --pos 5226924 \
      --ref-base C --alt-base G --name IVS_I_5 \
      --out-fa sim_ivs15.fa --out-vcf truth_ivs15.vcf --genotype 0/1
"""
import argparse, re, sys

def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None: recs.append((name, "".join(buf)))
            name, buf = line[1:], []
        else:
            buf.append(line.strip())
    if name is not None: recs.append((name, "".join(buf)))
    return recs

def write_fasta(path, recs, width=60):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(">" + name + "\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i+width] + "\n")

def parse_offset(header):
    m = re.match(r"^(.+):(\d+)-(\d+)$", header.split()[0])
    return (m.group(1), int(m.group(2))) if m else (header.split()[0], 1)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--pos", type=int, required=True, help="1-based, full-ref coords")
    p.add_argument("--ref-base", required=True)
    p.add_argument("--alt-base", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--out-vcf", required=True)
    p.add_argument("--genotype", default="0/1")
    a = p.parse_args()

    recs, out, done = read_fasta(a.ref), [], False
    for header, seq in recs:
        chrom, off = parse_offset(header)
        if chrom != a.chrom:
            out.append((header, seq)); continue
        idx = a.pos - off
        if not (0 <= idx < len(seq)):
            sys.exit("ERROR: pos %d outside contig %s" % (a.pos, header))
        actual = seq[idx].upper()
        if actual != a.ref_base.upper():
            sys.stderr.write("WARNING: reference base at %s:%d is %s, not %s "
                             "(proceeding, but check your coordinates)\n"
                             % (a.chrom, a.pos, actual, a.ref_base))
        new = seq[:idx] + a.alt_base + seq[idx+1:]
        out.append((header, new)); done = True
        sys.stderr.write("%s:%d  %s->%s  (contig index %d)\n"
                         % (a.chrom, a.pos, actual, a.alt_base, idx))
    if not done:
        sys.exit("ERROR: no contig matched %s" % a.chrom)

    write_fasta(a.out_fa, out)
    with open(a.out_vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttruth\n")
        fh.write("%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\tGT\t%s\n"
                 % (a.chrom, a.pos, a.name, a.ref_base, a.alt_base, a.genotype))
    sys.stderr.write("wrote %s and %s\n" % (a.out_fa, a.out_vcf))

if __name__ == "__main__":
    main()
    