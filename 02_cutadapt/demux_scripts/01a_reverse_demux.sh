#!/usr/bin/env bash --login
set -euo pipefail # Kills script in case of command failure

# Usage: 01a_reverse_demux.sh <SAMPLE>

# ==================== SETUP =========================

SAMPLE=$1

# Raw data
RAW_DIR=$(cat ../../00_prep/raw_dir.txt)

# Reverse primer fasta
PRIMER_REV=../primers/reverse_primers.fa

# Extract primer names from fasta file
mapfile -t PRIMERS < <(grep '^>' "$PRIMER_REV" | sed 's/^>//' | tr -d '\r')

# Append unknown
PRIMERS+=("unknown")

# Input directories
IN_R1=${RAW_DIR}/${SAMPLE}_R1_001.fastq.gz
IN_R2=${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz

# Output directories
OUTDIR=../demux_out/01a_reverse/${SAMPLE}
HCODIR=../demux_out/01a_reverse/${SAMPLE}/HCO2198
SSUDIR=../demux_out/01a_reverse/${SAMPLE}/SSUR22
UNKDIR=../demux_out/01a_reverse/${SAMPLE}/unknown
MANDIR=../demux_out/01a_reverse/summaries/manifests
LOGDIR=../demux_out/01a_reverse/summaries/CA_logs
SUMDIR=../demux_out/01a_reverse/summaries/summary_data

mkdir -p "$OUTDIR" "$MANDIR" "$LOGDIR" "$SUMDIR" \
"$HCODIR" "$SSUDIR" "$UNKDIR"

# ====================== CUTADAPT ===========================

echo "Processing reverse demux for $SAMPLE"

# Run cutadapt, bin by reverse primers
# Use --action=none if you want to skip the trimming of adapters and only demultiplex
cutadapt \
  -j 4 \
  -g file:"${PRIMER_REV}" \
  --pair-filter=any \
  --info-file="${LOGDIR}/${SAMPLE}_rev_info.tsv" \
  -o "${OUTDIR}/${SAMPLE}_{name}_R2.fastq.gz" \
  -p "${OUTDIR}/${SAMPLE}_{name}_R1.fastq.gz" \
  "${IN_R2}" "${IN_R1}" > $LOGDIR/${SAMPLE}_reverse.log 2>&1

### IMPORTANT NOTE!!!!!!
# During demultiplexing, cutadapt only uses forward primers to decide where each read is written.
# To 'cheat' this, we switch the positions of the reverse primers, reverse read input, and reverse read output
# So we 1) change -G to -g, 
# 2) switch -o and -p inputs (-o is technically for 'forward' reads).
# 3) switch the orders of the input files so the reverse reads come first.
# We do NOT do this in steps that use forward primers.

status=$? # Cutadapt exit status

# Move output files to corresponding rev marker subdirectory
mv ${OUTDIR}/*HCO*.* $HCODIR
mv ${OUTDIR}/*SSU*.* $SSUDIR
mv ${OUTDIR}/*unk*.* $UNKDIR


# ====================== CREATING MANIFEST AND SUMMARY FILES ===========================

echo "Creating manifest and summary files for $SAMPLE"

# Helper functions to count reads and calculate average read length
count_reads(){ echo $(( $(zcat < "$1" | wc -l) / 4 )); } # Opens file, counts lines, divides by 4

avg_len() {
    zcat < "$1" | awk 'NR % 4 == 2 { total += length($0); count++ }
                     END { print total / count }'
} # Opens file, captures every 2nd line in each 4-line chunk, 
# incrementally sums length of line, incrementally counts number of lines processed, calculates average


# Setting manifest and summary files
outMan="${MANDIR}/${SAMPLE}_rev_out_man.csv"
outSum="${SUMDIR}/${SAMPLE}_rev_out_sum.csv"

# Creating headers
man_header="sample,reverse_bin,out_R1_path,out_R2_path,CA_exit_status"
sum_header="sample,reverse_bin,readCount,FR_countsMatch,binSum_countsMatch,r1_avgLen,r2_avgLen"

printf "%s\n" "$man_header" > "$outMan"
printf "%s\n" "$sum_header" > "$outSum"

# Adding line for raw files
abs_raw_R1=$(realpath ${IN_R1}) # Getting raw paths
abs_raw_R2=$(realpath ${IN_R2})

raw_manline="${SAMPLE},raw,${abs_raw_R1},${abs_raw_R2},$status"
printf "%s\n" "$raw_manline" >> "$outMan"

# Raw data read counts
count_raw_r=$(count_reads "$IN_R2" 2>/dev/null || echo 0) # Getting raw counts and readLengths
count_raw_f=$(count_reads "$IN_R1" 2>/dev/null || echo 0)
avg_raw_f=$(avg_len "$IN_R1" 2>/dev/null || echo 0)
avg_raw_r=$(avg_len "$IN_R2" 2>/dev/null || echo 0)

if [[ $count_raw_f != $count_raw_r ]]; then
  FR_countsMatch="NO"
else
  FR_countsMatch="YES"
fi

raw_sumLine="${SAMPLE},raw,${count_raw_r},${FR_countsMatch},,${avg_raw_f},${avg_raw_r}"
printf "%s\n" "$raw_sumLine" >> "$outSum"

# Initalize count sums.
count_rev_total_r=0
count_rev_total_f=0

# Build each line by adding values to each column in primer order
for primer in "${PRIMERS[@]}"; do       # For each primer:

  # Get paths and output to manifest file
  r1_out_path=$(realpath ${OUTDIR}/${primer}/*R1.*)
  r2_out_path=$(realpath ${OUTDIR}/${primer}/*R2.*)
  manLine="${SAMPLE},${primer},${r1_out_path},${r2_out_path},${status}"
  printf "%s\n" "$manLine" >> "$outMan"

  # Get counts and average lengths
  count_f=$(count_reads "$r1_out_path" 2>/dev/null || echo 0)
  count_r=$(count_reads "$r2_out_path" 2>/dev/null || echo 0)
  avg_f=$(avg_len "$r1_out_path" 2>/dev/null || echo 0)
  avg_r=$(avg_len "$r2_out_path" 2>/dev/null || echo 0)

  # Running total of count sums
  count_rev_total_r=$((count_rev_total_r + count_r))
  count_rev_total_f=$((count_rev_total_f + count_f))

  # Check forward/reverse per-primer; if any mismatch, flag as NO
  if [[ $count_r != $count_f ]]; then
    FR_countsMatch="NO"
  else
    FR_countsMatch="YES"
  fi

  # Build line and output to summary file
  sumLine="${SAMPLE},${primer},${count_r},${FR_countsMatch},,${avg_f},${avg_r}"
  printf "%s\n" "$sumLine" >> "$outSum"
done

# Check if read totals match
if [[ $count_rev_total_r == $count_raw_r && $count_rev_total_f == $count_raw_f ]]; then
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
# Keep header unchanged, ensure at least 5 fields exist, only fill if empty, fill with Y/N check
