#!/bin/bash
# Phasing accuracy vs HPRC assembly-derived truth (whatshap 2.8).
# Truth VCFs carry sample name 'syndip', pipeline VCFs 'SAMPLE'; both are
# reheadered to the sample ID so whatshap can pair them.
set -e
module load BCFtools
WHATSHAP=~/whatshap_env/bin/whatshap
mkdir -p phasing
for S in HG02071 HG02083 HG02514 HG02074 HG02622; do
  echo "$S" > /tmp/sn.txt
  bcftools reheader -s /tmp/sn.txt wide_val/$S/truth.norm.vcf.gz -o phasing/$S.truth.vcf.gz
  bcftools index -t phasing/$S.truth.vcf.gz
  bcftools reheader -s /tmp/sn.txt wide_val/$S/pipeline.norm.vcf.gz -o phasing/$S.pipeline.vcf.gz
  bcftools index -t phasing/$S.pipeline.vcf.gz
  $WHATSHAP compare --names truth,pipeline \
    --tsv-pairwise phasing/$S.compare.tsv \
    phasing/$S.truth.vcf.gz phasing/$S.pipeline.vcf.gz
done
