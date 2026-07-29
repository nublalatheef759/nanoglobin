#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --output=pipe_%x_%j.log

S=$1
source ~/.bashrc 2>/dev/null
conda activate snakemake
cd /rds/projects/e/elhamsak-nanothal

sed "s/HG02071/${S}/g" config_hprc.yml > config_${S}.yml

snakemake --unlock --configfile config_${S}.yml 2>/dev/null
snakemake --configfile config_${S}.yml --use-conda --cores 16 \
  --set-threads minimap_align=16 -p \
  results/${S}.globin_only.vcf

echo "=== $S PIPELINE DONE $(date) ==="
ls -lh results/${S}.globin_only.vcf 2>/dev/null || echo "VCF NOT PRODUCED — check log"
