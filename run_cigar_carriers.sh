#!/bin/bash
# Run with:  source env_sv.sh && ./run_cigar_carriers.sh
# env_sv.sh loads bear-apps/2022a, which sets the library paths samtools
# needs. Calling the samtools binary by absolute path without the module
# stack fails on libbz2 and silently returns zero reads for every sample.
# CIGAR deletion typing on DRAGEN-labelled HPRC carriers and normal controls.
# Genotype labels are DRAGEN comparator calls, not independent truth.
set -e
SAMTOOLS=$(which samtools)
mkdir -p evaluation
for entry in "NA21106:-a4.2/aa" "HG00642:-a3.7/aa" "HG03136:-a3.7/-a3.7" \
             "HG03862:--/aa" "HG00735:aaa3.7/aa" "HG02071:aa/aa" "HG02074:aa/aa"; do
  s="${entry%%:*}"; gt="${entry##*:}"
  echo "=== $s ($gt) ==="
  python scripts/detect_deletions_cigar.py --bam sorted_reads/$s.bam \
    --name "$s" --genotype="$gt" --samtools "$SAMTOOLS" --cnvs databases/cnvs.csv
  echo ""
done
