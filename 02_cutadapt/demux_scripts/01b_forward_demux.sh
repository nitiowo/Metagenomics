#!/usr/bin/env bash --login
set -euo pipefail

# Usage: forward_demux.sh <SAMPLE> <REV_BIN>

# Forward primer fasta
PRIMER_FWD=../primers/fwd_primers.fa

# Extract primer names from fasta file
mapfile -t PRIMERS < <(grep '^>' "$PRIMER_FWD" | sed 's/^>//' | tr -d '\r')

# Append unknown
PRIMERS+=("unknown")

# Define paths
SAMPLE=$1
REV_BIN=$2
COMBO="${SAMPLE}_${REV_BIN}"

# Specific input folder for each sample x reverse primer combination:
IN_SamPrim_DIR=../demux_out/01a_reverse/${SAMPLE}/${REV_BIN} 

# Input fastq.gz files
IN_R1=${IN_SamPrim_DIR}/${COMBO}_R1.fastq.gz
IN_R2=${IN_SamPrim_DIR}/${COMBO}_R2.fastq.gz

# Output directories
OUT_FWD_DIR=../demux_out/01b_forward/${SAMPLE}/${REV_BIN}
MANDIR=../demux_out/01b_forward/summaries/manifests
LOGDIR=../demux_out/01b_forward/summaries/CA_logs
SUMDIR=../demux_out/01b_forward/summaries/summary_data

# Separate output directories for each fwd primer (and unknown)
LCODIR=${OUT_FWD_DIR}/LCO1490
MCODIR=${OUT_FWD_DIR}/mICOintF
SSUDIR=${OUT_FWD_DIR}/SSUF04
UNKDIR=${OUT_FWD_DIR}/unknown

mkdir -p "$OUT_FWD_DIR" "$SUMDIR" "$LOGDIR" "$MANDIR" \
"$LCODIR" "$MCODIR" "$SSUDIR" "$UNKDIR"

echo "Processing forward demux for $SAMPLE"

# Activate cutadapt through conda or module if necessary
eval "$(conda shell.bash hook)"
conda activate cutadaptenv

# Run cutadapt, bin by reverse primers
# Use --action=none if you want to skip the trimming of adapters and only demultiplex
cutadapt \
  -j 4 \
  -g file:"${PRIMER_FWD}" \
  --pair-filter=any \
  --info-file="${LOGDIR}/${COMBO}_fwd_info.tsv" \
  -o "${OUT_FWD_DIR}/${COMBO}_{name}_R1.fastq.gz" \
  -p "${OUT_FWD_DIR}/${COMBO}_{name}_R2.fastq.gz" \
  "${IN_R1}" "${IN_R2}" > $LOGDIR/${COMBO}_fwd.log 2>&1

status=$? # Cutadapt exit status

# Move output files to corresponding fwd marker subdirectory
mv ${OUT_FWD_DIR}/*LCO*.* $LCODIR
mv ${OUT_FWD_DIR}/*mICO*.* $MCODIR
mv ${OUT_FWD_DIR}/*SSUR*.* $SSUDIR
mv ${OUT_FWD_DIR}/*unk*.* $UNKDIR


# ====================== CREATING MANIFEST AND SUMMARY FILES ===========================

echo "Creating manifest and summary files for $SAMPLE"

# Helper functions to count reads and calculate average read length
count_reads(){ echo $(( $(zcat < "$1" | wc -l) / 4 )); } # Opens file, counts lines, divides by 4

avg_len() {
    zcat < "$1" | awk 'NR % 4 == 2 { total += length($0); count++ }
                     END { print total / count }'
} # Opens file, captures every 2nd line in each 4-line chunk, 
# Incrementally sums length of line, incrementally counts number of lines processed, calculates average

# Get raw files
RAW_DIR=$(cat ../../00_prep/raw_dir.txt)

Raw_R1=${RAW_DIR}/${SAMPLE}_R1_001.fastq.gz
Raw_R2=${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz

# Setting manifest and summary files
outMan="${MANDIR}/${COMBO}_fwd_out_man.csv"
outSum="${SUMDIR}/${COMBO}_fwd_out_sum.csv"

# Manifest header
samp_man="../demux_out/01a_reverse/summaries/manifests/${SAMPLE}_rev_out_man.csv"
samp_man_head=$(head -n 1 "$samp_man") 
samp_man_head+=",fwd_bin"
printf "%s\n" "$samp_man_head" > "$outMan"

# Manifest raw file line
samp_man_raw=$(head -n 2 "$samp_man" | tail -n 1)
samp_man_raw+=",raw"
printf "%s\n" "$samp_man_raw" >> "$outMan"

# Manifest parent revbin line
samp_man_revbin=$(grep ${REV_BIN} ${samp_man})
samp_man_revbin+=",parent"
printf "%s\n" "$samp_man_revbin" >> "$outMan"

# Summary file header
samp_sum="../demux_out/01a_reverse/summaries/summary_data/${SAMPLE}_rev_out_sum.csv"
samp_sum_head=$(head -n 1 "$samp_sum") 
samp_sum_head+=",fwd_bin"
printf "%s\n" "$samp_sum_head" > "$outSum"

# Summary raw file line
samp_sum_raw=$(head -n 2 "$samp_sum" | tail -n 1)
samp_sum_raw+=",raw"
printf "%s\n" "$samp_sum_raw" >> "$outSum"

# Summary parent revbin line
samp_sum_revbin=$(grep ${REV_BIN} ${samp_sum})
samp_sum_revbin+=",parent"
printf "%s\n" "$samp_sum_revbin" >> "$outSum"

# Get parent revbin count data
parentBin_count=$(echo ${samp_sum_revbin} | cut -d ',' -f 3)

# Initalize count sums.
count_fwd_total_r=0
count_fwd_total_f=0

# Build each line by adding values to each column in primer order
for primer in "${PRIMERS[@]}"; do       # For each primer:

  # Get paths, and output to manifest file
  r1_out_path=$(realpath ${OUT_FWD_DIR}/${primer}/*R1.*)
  r2_out_path=$(realpath ${OUT_FWD_DIR}/${primer}/*R2.*)
  manLine="${SAMPLE},${REV_BIN},${r1_out_path},${r2_out_path},${status},${primer}"
  printf "%s\n" "$manLine" >> "$outMan"

  # Get counts and average lengths
  count_f=$(count_reads "$r1_out_path" 2>/dev/null || echo 0)
  count_r=$(count_reads "$r2_out_path" 2>/dev/null || echo 0)
  avg_f=$(avg_len "$r1_out_path" 2>/dev/null || echo 0)
  avg_r=$(avg_len "$r2_out_path" 2>/dev/null || echo 0)

  # Running total of count sums
  count_fwd_total_r=$((count_fwd_total_r + count_r))
  count_fwd_total_f=$((count_fwd_total_f + count_f))

  # Check forward/reverse per-primer; if any mismatch, flag as NO
  if [[ $count_r != $count_f ]]; then
    FR_countsMatch="NO"
  else
    FR_countsMatch="YES"
  fi

  # Build line and output to summary file
  sumLine="${SAMPLE},${REV_BIN},${count_r},${FR_countsMatch},,${avg_f},${avg_r},${primer}"
  printf "%s\n" "$sumLine" >> "$outSum"
done

# Check if read totals match
if [[ $count_fwd_total_r == $count_fwd_total_f && $count_fwd_total_f == $parentBin_count ]]; then
  Tot_countsMatch="YES"
else
  Tot_countsMatch="NO"
fi

# Add total counts match Y/N to every line in summary file
tmp=$(mktemp "${outSum}.tmp.XXXX")
awk -F, -v OFS=, -v val="$Tot_countsMatch" '
  NR==1 { print; next }                  
  {
    for(i=NF+1; i<=5; i++) $i = ""
    if ($5 == "") $5 = val               
    print
  }
' "$outSum" > "$tmp" && mv "$tmp" "$outSum"
# Keep header unchanged, ensure at least 5 fields exist, only fill 5th field if empty, fill with Y/N check