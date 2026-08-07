#!/bin/bash
# Alpha-cluster gene-body tier: MANE Select spans (HBZ, HBA2, HBA1).
# Result: 5 SNV TP / 0 FN / 0 FP; no indels assessed. Small denominator --
# diagnostic only, confirms indel misses lie outside globin gene bodies.
set -e
for S in HG02071 HG02083 HG02514 HG02074 HG02622; do
  singularity exec happy.sif /opt/hap.py/bin/hap.py \
    wide_val/$S/truth.norm.vcf.gz wide_val/$S/pipeline.norm.vcf.gz \
    -r reference/hg38.fa -T regions/alpha_gene_bodies.bed --threads 4 \
    -o wide_val/$S/happy_alpha
done
