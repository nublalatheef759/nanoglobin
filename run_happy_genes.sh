#!/bin/bash
# Beta-cluster gene-body tier: MANE Select spans (HBB, HBD, HBG1, HBG2, HBE1).
# Result: 132 SNV TP / 1 FN / 1 FP; 21 INDEL TP / 0 FN / 0 FP across 5 samples.
set -e
for S in HG02071 HG02083 HG02514 HG02074 HG02622; do
  singularity exec happy.sif /opt/hap.py/bin/hap.py \
    wide_val/$S/truth.norm.vcf.gz wide_val/$S/pipeline.norm.vcf.gz \
    -r reference/hg38.fa -T regions/beta_gene_bodies.bed --threads 4 \
    -o wide_val/$S/happy_genes_v2
done
