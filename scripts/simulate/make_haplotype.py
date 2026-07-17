#!/usr/bin/env python3
"""
make_haplotype.py -- apply a deletion to a simulation reference and emit a truth VCF.

The simulation reference is assumed to have been cut with samtools faidx, so its
contigs are named like 'chr16:165000-182000'. Coordinates on the command line are
given in FULL reference space (e.g. chr16:173200) and translated internally.

Usage:
  python make_haplotype.py --ref sim_target.fa --chrom chr16 \
      --del-start 173200 --del-end 176899 --name a37 \
      --out-fa sim_a37.fa --out-vcf truth_a37.vcf --genotype 0/1
"""
import argparse
import re
import sys


def read_fasta(path):
    """Return list of (header, sequence) preserving order."""
    recs = []
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path, recs, width=60):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(">" + name + "\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def parse_contig(header):
    """'chr16:165000-182000' -> ('chr16', 165000, 182000). Bare 'chr16' -> ('chr16', 1, None)."""
    tok = header.split()[0]
    m = re.match(r"^(.+):(\d+)-(\d+)$", tok)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return tok, 1, None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True)
    p.add_argument("--chrom", required=True, help="chromosome in FULL reference space, e.g. chr16")
    p.add_argument("--del-start", type=int, required=True, help="1-based inclusive, full-ref coords")
    p.add_argument("--del-end", type=int, required=True, help="1-based inclusive, full-ref coords")
    p.add_argument("--name", required=True, help="variant label, e.g. a37")
    p.add_argument("--out-fa", required=True)
    p.add_argument("--out-vcf", required=True)
    p.add_argument("--genotype", default="1/1")
    a = p.parse_args()

    if a.del_end < a.del_start:
        sys.exit("ERROR: --del-end is before --del-start")

    recs = read_fasta(a.ref)
    if not recs:
        sys.exit("ERROR: no sequences read from " + a.ref)

    out, hits = [], 0
    anchor_base = None
    for header, seq in recs:
        chrom, off, end = parse_contig(header)
        if chrom != a.chrom:
            out.append((header, seq))
            continue

        # contig position 1 == full-ref position `off`, so 0-based index = P - off
        lo = a.del_start - off
        hi = a.del_end - off          # inclusive
        if lo < 1 or hi >= len(seq):
            sys.exit(
                "ERROR: deletion %d-%d falls outside contig %s (covers %d-%d)"
                % (a.del_start, a.del_end, header, off, off + len(seq) - 1)
            )

        anchor_base = seq[lo - 1]      # base immediately before the deletion
        new = seq[:lo] + seq[hi + 1:]
        out.append((header, new))
        hits += 1
        sys.stderr.write(
            "%s: %d bp -> %d bp (removed %d bp at %s:%d-%d)\n"
            % (header, len(seq), len(new), len(seq) - len(new),
               a.chrom, a.del_start, a.del_end)
        )

    if hits == 0:
        sys.exit("ERROR: no contig matched chromosome '%s'" % a.chrom)
    if hits > 1:
        sys.exit("ERROR: %d contigs matched '%s'; ambiguous" % (hits, a.chrom))

    write_fasta(a.out_fa, out)

    svlen = -(a.del_end - a.del_start + 1)
    pos = a.del_start - 1              # VCF anchors on the base before the event
    with open(a.out_vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        fh.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT">\n')
        fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
        fh.write('##ALT=<ID=DEL,Description="Deletion">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttruth\n")
        fh.write(
            "%s\t%d\t%s\t%s\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=%d;SVLEN=%d\tGT\t%s\n"
            % (a.chrom, pos, a.name, anchor_base, a.del_end, svlen, a.genotype)
        )

    sys.stderr.write("wrote %s and %s\n" % (a.out_fa, a.out_vcf))


if __name__ == "__main__":
    main()
