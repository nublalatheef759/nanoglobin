#!/bin/bash
# run_happy.sh -- GA4GH-standard concordance with hap.py, via Singularity
# (sidesteps the conda pysam/libdeflate conflict entirely).
#
# Runs HG002 first (GIAB is exactly what hap.py is built for), then the 5 HPRC
# samples, comparing the pipeline's globin calls against the same truth VCFs used
# for the bcftools-isec validation. Prints hap.py's summary so you can compare
# its TP/FN/FP to your isec numbers.
#
# Usage:  bash run_happy.sh
# (login node for the pull; the hap.py runs are quick on the globin region)

set -e
PROJ=/rds/projects/e/elhamsak-nanothal
cd $PROJ
REF=$PROJ/reference/hg38.fa
SIF=$PROJ/happy.sif

module purge
module load bb-singularity-conf/live 2>/dev/null || true

# --- container already pulled as happy.sif (jmcdani20/hap.py:v0.3.12) ---
if [ ! -s "$SIF" ]; then
  echo "=== pulling hap.py container $(date) ==="
  singularity pull "$SIF" docker://jmcdani20/hap.py:v0.3.12
fi

# hap.py inside the jmcdani20 image lives at /opt/hap.py/bin/hap.py
HAPPY=/opt/hap.py/bin/hap.py

run_happy () {
  local S=$1 TRUTH=$2 QUERY=$3 BED=$4
  echo ""
  echo "=== hap.py: $S $(date) ==="
  mkdir -p happy_out/$S
  local bedarg=""
  [ -n "$BED" ] && [ -s "$BED" ] && bedarg="-T $BED"
  singularity exec "$SIF" $HAPPY \
    "$TRUTH" "$QUERY" \
    -r "$REF" \
    -o happy_out/$S/happy \
    $bedarg 2>&1 | tail -6
  echo "--- $S summary ---"
  [ -s happy_out/$S/happy.summary.csv ] && \
    column -t -s, happy_out/$S/happy.summary.csv | head -6
}

# --- HG002 first (cleanest test; GIAB high-conf globin) ---
# uses the normalised truth + pipeline VCFs already built in hprc_validation/HG002
if [ -s hprc_validation/HG002/truth_globin.norm.vcf.gz ]; then
  run_happy HG002 \
    hprc_validation/HG002/truth_globin.norm.vcf.gz \
    hprc_validation/HG002/pipeline_globin.norm.vcf.gz \
    hprc_validation/HG002/globin_highconf.bed
fi

# --- the 5 HPRC SNV samples ---
for S in HG02071 HG02083 HG02514 HG02074 HG02622; do
  T=hprc_validation/$S/truth.norm.vcf.gz
  Q=hprc_validation/$S/pipeline.norm.vcf.gz
  if [ -s "$T" ] && [ -s "$Q" ]; then
    run_happy $S "$T" "$Q" ""
  else
    echo "skip $S (normalised VCFs not found at $T / $Q)"
  fi
done

echo ""
echo "=== DONE. Compare hap.py TP/FN/FP below to your isec numbers ==="
echo "Per-sample summaries are in happy_out/<sample>/happy.summary.csv"
for S in HG002 HG02071 HG02083 HG02514 HG02074 HG02622; do
  f=happy_out/$S/happy.summary.csv
  if [ -s "$f" ]; then
    echo "--- $S (INDEL + SNP rows) ---"
    column -t -s, "$f" | awk 'NR==1 || /PASS/'
  fi
done
