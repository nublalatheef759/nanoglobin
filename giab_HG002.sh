#!/bin/bash
#SBATCH --job-name=giab_HG002
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=giab_HG002_%j.log

# GIAB HG002 gold-standard SNV/indel validation.
# HG002's HBA is normal (aa/aa) and GIAB v4.2.1 excludes segmental duplications,
# so this primarily validates HBB (chr11) against the community-standard truth.

S=HG002
source ~/.bashrc 2>/dev/null
PROJ=/rds/projects/e/elhamsak-nanothal
cd $PROJ
mkdir -p scratch_tmp hprc_validation/$S
export TMPDIR=$PROJ/scratch_tmp

CRAM="s3://ont-open-data/giab_2023.05/analysis/hg002/sup/PAO83395.pass.cram"

# ---------- 1. extract chr11 + chr16 reads from the aligned CRAM ----------
if [ ! -s fastq/${S}.fastq ]; then
  echo "=== $S: extracting chr11+chr16 from GIAB CRAM $(date) ==="
  module purge; module load bear-apps/2022a/live; module load SAMtools/1.16.1-GCC-11.3.0
  samtools view -T reference/hg38.fa -b "$CRAM" chr11 chr16 \
    > scratch_tmp/${S}.regions.bam
  echo "extracted reads: $(samtools view -c scratch_tmp/${S}.regions.bam)"
  # aligned BAM -> FASTQ (raw reads, so YOUR pipeline re-aligns = clean test)
  samtools fastq -@ 16 scratch_tmp/${S}.regions.bam > fastq/${S}.fastq
  echo "$S FASTQ reads: $(( $(wc -l < fastq/${S}.fastq)/4 ))"
  rm -f scratch_tmp/${S}.regions.bam
else
  echo "=== $S: FASTQ exists, skipping extraction ==="
fi

# ---------- 2. run the pipeline ----------
echo "=== $S: pipeline $(date) ==="
module purge; source ~/.bashrc 2>/dev/null; conda activate snakemake
cd $PROJ
export TMPDIR=$PROJ/scratch_tmp
sed "s/HG02071/${S}/g" config_hprc.yml > config_${S}.yml
snakemake --unlock --configfile config_${S}.yml 2>/dev/null
snakemake --configfile config_${S}.yml --use-conda --cores 16 \
  --set-threads minimap_align=16 -p \
  results/${S}.globin_only.vcf 2>&1
echo "pipeline calls by chr:"
grep -v "^#" results/${S}.globin_only.vcf | cut -f1 | sort | uniq -c

# ---------- 3. download NIST v4.2.1 truth (VCF + high-confidence BED) ----------
echo "=== $S: downloading NIST truth $(date) ==="
cd $PROJ/hprc_validation/$S
BASE="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
[ ! -s truth.vcf.gz ] && curl -sL "$BASE/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" -o truth.vcf.gz
[ ! -s highconf.bed ] && curl -sL "$BASE/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" -o highconf.bed
echo "truth.vcf.gz: $(ls -lh truth.vcf.gz | awk '{print $5}')"
echo "highconf.bed: $(ls -lh highconf.bed | awk '{print $5}')"

# ---------- 4. compare in globin regions INTERSECTED with high-conf BED ----------
echo "=== $S: comparison $(date) ==="
module purge; module load bear-apps/2022a/live
module load BCFtools/1.15.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.3.0 2>/dev/null

# globin regions of interest
printf "chr11\t5225000\t5228000\nchr16\t170000\t178000\n" > globin.bed

# high-confidence globin = globin regions ∩ GIAB high-conf BED
bedtools intersect -a globin.bed -b highconf.bed > globin_highconf.bed
echo "high-confidence globin regions:"
cat globin_highconf.bed
echo "bp of high-conf globin: $(awk '{s+=$3-$2} END{print s}' globin_highconf.bed)"

# index + normalise truth, restrict to high-conf globin
bcftools index -t truth.vcf.gz 2>/dev/null
bcftools view -R globin_highconf.bed truth.vcf.gz 2>/dev/null \
  | bcftools norm -f $PROJ/reference/hg38.fa 2>/dev/null \
  | bcftools view -Oz -o truth_globin.norm.vcf.gz
bcftools index -t truth_globin.norm.vcf.gz

# normalise pipeline calls, restrict to same high-conf globin
bcftools norm -f $PROJ/reference/hg38.fa $PROJ/results/${S}.globin_only.vcf 2>/dev/null \
  | bcftools view -R globin_highconf.bed -Oz -o pipeline_globin.norm.vcf.gz 2>/dev/null
bcftools index -t pipeline_globin.norm.vcf.gz

# intersect
rm -rf isec; bcftools isec truth_globin.norm.vcf.gz pipeline_globin.norm.vcf.gz -p isec
TP=$(grep -vc "^#" isec/0002.vcf)
FN=$(grep -vc "^#" isec/0000.vcf)
FP=$(grep -vc "^#" isec/0001.vcf)

echo "============================================="
echo "GIAB HG002 (high-confidence globin regions):"
echo "  TP=$TP  FN=$FN  FP=$FP"
[ $((TP+FN)) -gt 0 ] && echo "  Sensitivity: $(echo "scale=1;$TP*100/($TP+$FN)"|bc)%"
[ $((TP+FP)) -gt 0 ] && echo "  Precision:   $(echo "scale=1;$TP*100/($TP+$FP)"|bc)%"
echo "============================================="
echo "=== DONE $(date) ==="
