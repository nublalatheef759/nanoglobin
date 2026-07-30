#!/bin/bash
#SBATCH --job-name=dipcall_carriers
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=dipcall_carriers_%j.log
#
# dipcall_carriers.sh -- run the existing dipcall workflow on the three alpha-
# deletion carriers, to obtain assembly-based truth for the alpha-globin region
# (independent of DRAGEN). Reuses dipcall_one.sh's exact steps.
#
# The alpha deletions live in a segmental duplication (HBA1/HBA2), so the key
# question is whether the Verkko hybrid assemblies resolve the deletion cleanly.
# The per-sample GLOBIN TRUTH block at the end prints what the assembly shows.
#
# Run:  sbatch dipcall_carriers.sh
# (downloads are large; the job downloads then dipcalls each carrier in turn)

BASE_S3="s3://human-pangenomics/submissions/bed38f3f-cac4-4703-9852-181bc362cea1--R3_verkko-v2.3.2_hybrid_assembly"
REF=/rds/projects/e/elhamsak-nanothal/reference/hg38.fa
ROOT=/rds/projects/e/elhamsak-nanothal/hprc_validation
MM=/rds/homes/f/fnl759/.conda/envs/dipcall/bin/minimap2
K8=/rds/homes/f/fnl759/.conda/envs/dipcall/bin/k8
PAFTOOLS=/rds/homes/f/fnl759/.conda/envs/dipcall/bin/paftools.js
BEDTK=/rds/homes/f/fnl759/.conda/envs/dipcall/bin/bedtk

# the three deletion carriers (DRAGEN-genotyped: HG00642 -a3.7/aa,
# HG03136 -a3.7/-a3.7, NA21106 -a4.2/aa)
for S in HG00642 HG03136 NA21106; do
  echo ""
  echo "############ $S  $(date) ############"
  DIR=$ROOT/$S
  mkdir -p $DIR
  cd $DIR

  # --- download both haplotype assemblies (skip if present) ---
  if [ ! -s ${S}.assembly.refOriented.haplotype1.fasta ] || \
     [ ! -s ${S}.assembly.refOriented.haplotype2.fasta ]; then
    ( module purge
      module load bear-apps/2023a/live
      module load awscli/2.17.54-GCCcore-12.3.0
      for H in 1 2; do
        echo "=== $S: downloading haplotype $H $(date) ==="
        aws s3 cp --no-sign-request \
          "$BASE_S3/$S/verkko-hi-c/${S}.assembly.refOriented.haplotype${H}.fasta" \
          ${S}.assembly.refOriented.haplotype${H}.fasta
      done )
  else
    echo "=== $S: assemblies already present, skipping download ==="
  fi

  # if download failed (sample may use a different assembler subdir), note + skip
  if [ ! -s ${S}.assembly.refOriented.haplotype1.fasta ]; then
    echo "!!! $S: assembly not found under verkko-hi-c/ -- check available paths:"
    ( module purge; module load bear-apps/2023a/live
      module load awscli/2.17.54-GCCcore-12.3.0
      aws s3 ls --no-sign-request "$BASE_S3/$S/" )
    continue
  fi

  # --- dipcall (exact reuse of dipcall_one.sh) ---
  source ~/.bashrc 2>/dev/null
  conda activate dipcall
  for H in 1 2; do
    $MM -c --paf-no-hit -xasm5 --cs -t16 $REF \
      ${S}.assembly.refOriented.haplotype${H}.fasta 2> ${S}.hap${H}.paf.gz.log \
      | gzip > ${S}.hap${H}.paf.gz
    $MM -a -xasm5 --cs -t16 $REF \
      ${S}.assembly.refOriented.haplotype${H}.fasta 2> ${S}.hap${H}.sam.gz.log \
      | gzip > ${S}.hap${H}.sam.gz
    gzip -dc ${S}.hap${H}.paf.gz | sort -k6,6 -k8,8n \
      | $K8 $PAFTOOLS call - 2> ${S}.hap${H}.var.gz.vst | gzip > ${S}.hap${H}.var.gz
    gzip -dc ${S}.hap${H}.var.gz | grep ^R | cut -f2- > ${S}.hap${H}.bed
    $K8 /rds/homes/f/fnl759/.conda/envs/dipcall/bin/dipcall-aux.js samflt ${S}.hap${H}.sam.gz \
      | samtools sort -m4G --threads 4 -o ${S}.hap${H}.bam -
  done
  $BEDTK isec -m ${S}.hap1.bed ${S}.hap2.bed > ${S}.dip.bed
  /rds/homes/f/fnl759/.conda/envs/dipcall/bin/htsbox pileup -q5 -evcf $REF \
    ${S}.hap1.bam ${S}.hap2.bam \
    | /rds/homes/f/fnl759/.conda/envs/dipcall/bin/htsbox bgzip > ${S}.pair.vcf.gz
  $K8 /rds/homes/f/fnl759/.conda/envs/dipcall/bin/dipcall-aux.js vcfpair ${S}.pair.vcf.gz \
    | /rds/homes/f/fnl759/.conda/envs/dipcall/bin/htsbox bgzip > ${S}.dip.vcf.gz

  # --- the key output: what does the assembly show at the alpha locus? ---
  echo "=== $S ALPHA-GLOBIN TRUTH (chr16:170000-178000) ==="
  echo -n "HBA variant records: "
  zcat ${S}.dip.vcf.gz | grep -v "^#" \
    | awk '$1=="chr16" && $2>=170000 && $2<=178000' | wc -l
  echo "-- variants + the dip.bed coverage over the alpha region --"
  zcat ${S}.dip.vcf.gz | grep -v "^#" \
    | awk '$1=="chr16" && $2>=170000 && $2<=178000 {print $1,$2,$4,$5,$7,$10}'
  echo "-- callable (dip.bed) intervals over chr16:170000-178000 --"
  awk '$1=="chr16" && $3>=170000 && $2<=178000' ${S}.dip.bed
  echo "############ $S DIPCALL DONE $(date) ############"
done

echo ""
echo "=== INTERPRETATION ==="
echo "For a -a3.7/-a4.2 deletion, expect either explicit DEL variant records in the"
echo "alpha region, OR a GAP in the dip.bed callable intervals where the deleted"
echo "segment is (the assembly haplotype simply lacks that sequence). A clean gap"
echo "matching the known breakpoints = assembly-confirmed deletion truth,"
echo "independent of DRAGEN."
