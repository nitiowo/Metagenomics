#!/usr/bin/env bash --login
set -euo pipefail

# Usage: forward_demux.sh <SAMPLE>

# User inputs:
adapter1="GGTCAACAAATCATAAAGATATTGG"  #LCO1490
adapter2="GGWACWGGWTGAACWGTWTAYCCYCC"  #mICOintF
adapter3="GCTTGTCTCAAAGATTAAGCC"  #SSUF04

# Define paths
SAMPLE=$1
REV_PRIMER=$2
COMBO="${SAMPLE}_${REV_PRIMER}"

# Specific input folder for each sample x reverse primer combination:
IN_SamPrim_DIR=../demux_out/01a_reverse/${SAMPLE}/${REV_PRIMER} 

# Input fastq.gz files
IN_R1=${IN_SamPrim_DIR}/${SAMPLE}_${REV_PRIMER}_R1.fastq.gz
IN_R1=${IN_SamPrim_DIR}/${SAMPLE}_${REV_PRIMER}_R1.fastq.gz

# Output directories
OUT_REV_DIR=../demux_out/01b_forward/${SAMPLE}/${REV_PRIMER}
SUMDIR=../demux_out/summaries/reads_n_counts
LOGDIR=../demux_out/summaries/fwd_logs/${SAMPLE}

# Separate output directories for each fwd primer (and unknown)
LCODIR=${OUT_REV_DIR}/LCO1490
MCODIR=${OUT_REV_DIR}/mICOintF
SSUDIR=${OUT_REV_DIR}/SSUF04
UNKDIR=${OUT_REV_DIR}/unknown

mkdir -p "$OUT_REV_DIR" "$SUMDIR" "$LOGDIR" \
"$LCODIR" "$MCODIR" "$SSUDIR"

# Helper functions to count reads and calculate average read length
count_reads(){ echo $(( $(zcat < "$1" | wc -l) / 4 )); }
avg_len() {
    zcat < "$1" | awk 'NR % 4 == 2 { total += length($0); count++ }
                     END { print total / count }'
}

echo "Processing forward demux for $SAMPLE"

# Activate cutadapt through conda or module if necessary
eval "$(conda shell.bash hook)"
conda activate cutadaptenv

# Run cutadapt, bin by reverse primers
# Use --action=none if you want to skip the trimming of adapters and only demultiplex
cutadapt \
  -j 4 \
  -g LCO1490=$adapter1 -g mICOintF=$adapter2 -g SSUF04=$adapter3 \
  --pair-filter=any \
  --info-file="${LOGDIR}/${COMBO}_fwd_info.tsv" \
  -o "${OUT_REV_DIR}/${COMBO}_{name}_R1.fastq.gz" \
  -p "${OUT_REV_DIR}/${COMBO}_{name}_R2.fastq.gz" \
  "${IN_R1}" "${IN_R2}" > $LOGDIR/${COMBO}_fwd.log 2>&1

# Move output files to corresponding fwd marker subdirectory
mv ${OUT_REV_DIR}/*LCO*.* $LCODIR
mv ${OUT_REV_DIR}/*mICO*.* $MCODIR
mv ${OUT_REV_DIR}/*SSU*.* $SSUDIR
mv ${OUT_REV_DIR}/*unk*.* $UNKDIR


# _________________ Summarizing outputs _____________________
















































