#!/usr/bin/env bash --login
set -euo pipefail # Kills script in case of command failure

# Usage: 01a_reverse_demux.sh <SAMPLE>

# ==================== SETUP =========================

SAMPLE=$1

# Reverse primer fasta
PRIMER_REV=../primers/reverse_primers.fa

# Extract primer names from fasta file
mapfile -t PRIMERS < <(grep '^>' "$PRIMER_REV" | sed 's/^>//' | tr -d '\r')

# Append unknown
PRIMERS+=("unknown")

# Input directories
IN_R1=../../raw_data/${SAMPLE}_R1_001.fastq.gz
IN_R2=../../raw_data/${SAMPLE}_R2_001.fastq.gz

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

# Activate cutadapt through conda or module if necessary
eval "$(conda shell.bash hook)"
conda activate cutadaptenv


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

# Write output manifest lines to file
outMan="${MANDIR}/${SAMPLE}_rev_out_man.csv"

# Creating header
man_header="sample,reverse_bin,out_R1_path,out_R2_path,CA_exit_status"
printf "%s\n" "$man_header" > "$outMan"

# Adding line for raw files
abs_raw_R1=$(realpath ${IN_R1})
abs_raw_R2=$(realpath ${IN_R2})
rawline="${SAMPLE},raw,${abs_raw_R1},${abs_raw_R2},$status"
printf "%s\n" "$rawline" >> "$outMan"

# Build each line by adding values to each column in primer order
for primer in "${PRIMERS[@]}"; do       # For each primer:
  line="$SAMPLE"
  r1_out_path=$(realpath ${OUTDIR}/${primer}/*R1.*)
  r2_out_path=$(realpath ${OUTDIR}/${primer}/*R2.*)
  line+=",${primer},${r1_out_path},${r2_out_path},$status"
  printf "%s\n" "$line" >> "$outMan"
done

# Rscript reverse_summary.R $outMan $SUMDIR