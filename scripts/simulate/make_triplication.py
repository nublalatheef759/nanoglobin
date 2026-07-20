#!/usr/bin/env python3
"""
make_triplication.py -- simulate a copy-number GAIN (e.g. alpha-globin
triplication anti-3.7 / anti-4.2) by inserting a tandem duplicate of a segment
into the reference, and emit a truth VCF.

Companion to make_haplotype.py (deletions). A deletion removes a segment so reads
map back at reduced depth; a triplication inserts an extra copy so reads map back
at INCREASED depth over the affected span (~1.5x for het aaa/aa, ~2x for aaaa/aa).

The duplicated segment is inserted in tandem immediately after the original, so
reads simulated from this modified reference, when aligned to the NORMAL
reference, pile up to 1.5x over [dup_start, dup_end].

Usage:
  python make_triplication.py --ref sim_target.fa --chrom chr16 \
      --dup-start 173000 --dup-end 177000 --name anti_a37 \
      --out-fa sim_triple.fa --out-vcf truth_triple.vcf --genotype 0/1
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

def parse_contig(header):
    tok = header.split()[0]
    m = re.match(r"^(.+):(\d+)-(\d+)$", tok)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return tok, 1, None

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--dup-start", type=int, required=True, help="1-based, full-ref coords")
    p.add_argument("--dup-end", type=int, required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--out-vcf", required=True)
    p.add_argument("--genotype", default="0/1")
    a = p.parse_args()

    recs, out, done = read_fasta(a.ref), [], False
    for header, seq in recs:
        chrom, off, _ = parse_contig(header)
        if chrom != a.chrom:
            out.append((header, seq)); continue
        s = a.dup_start - off
        e = a.dup_end - off
        if not (0 <= s < e <= len(seq)):
            sys.exit("ERROR: dup region %d-%d outside contig %s (len %d, offset %d)"
                     % (a.dup_start, a.dup_end, header, len(seq), off))
        segment = seq[s:e]
        # insert a tandem copy immediately after the original segment
        newseq = seq[:e] + segment + seq[e:]
        out.append((header, newseq)); done = True
        sys.stderr.write("%s:%d-%d duplicated (%d bp inserted in tandem)\n"
                         % (a.chrom, a.dup_start, a.dup_end, len(segment)))
    if not done:
        sys.exit("ERROR: no contig matched %s" % a.chrom)

    write_fasta(a.out_fa, out)
    svlen = a.dup_end - a.dup_start
    with open(a.out_vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        fh.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttruth\n")
        fh.write("%s\t%d\t%s\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=%d;SVLEN=%d\tGT\t%s\n"
                 % (a.chrom, a.dup_start, a.name, a.dup_end, svlen, a.genotype))
    sys.stderr.write("wrote %s and %s\n" % (a.out_fa, a.out_vcf))

if __name__ == "__main__":
    main()
    