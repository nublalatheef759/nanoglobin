#!/bin/bash
#SBATCH --job-name=pipe_HG02071
#SBATCH --time=6:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --output=pipe_HG02071_%j.log

source ~/.bashrc
conda activate snakemake
cd /rds/projects/e/elhamsak-nanothal

snakemake --configfile config_hprc.yml --use-conda --cores 16 \
  --set-threads minimap_align=16 -p \
  results/HG02071.globin_only.vcf
echo "PIPELINE DONE"
