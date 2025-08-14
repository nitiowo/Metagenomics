#!/usr/bin/env bash --login
set -euo pipefail # Kills script in case of command failure

# Usage: 01a_reverse_demux.sh <SAMPLE>

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
MANDIR=../demux_out/summaries/manifests
LOGDIR=../demux_out/summaries/CA_logs

mkdir -p "$OUTDIR" "$MANDIR" "$LOGDIR" "$HCODIR" "$SSUDIR" "$UNKDIR"

echo "Processing reverse demux for $SAMPLE"

# Activate cutadapt through conda or module if necessary
eval "$(conda shell.bash hook)"
conda activate cutadaptenv

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

# Write output manifest lines to file
outMan="${MANDIR}/${SAMPLE}_rev_out_man.tsv"
man_header="sample"$'\t'"reverse_bin"$'\t'"file_R1"$'\t'"file_R2"$'\t'"CA_exit_status"    # Create header line
printf "%s\n" "$man_header" > "$outMan"

# Build each line by adding values to each column in primer order
line="$SAMPLE"                            # Value in first column
for primer in "${PRIMERS[@]}"; do       # For each primer:
  line="$SAMPLE"
  r1_out_path=$(realpath ${OUTDIR}/${primer}/*R1.*)
  r2_out_path=$(realpath ${OUTDIR}/${primer}/*R2.*)
  line+=$'\t'"${primer}"$'\t'"${r1_out_path}"$'\t'"${r2_out_path}"$'\t'"$status"  # Add count to the line (tab separated)
  printf "%s\n" "$line" >> "$outMan"
done

# Rscript reverse_summary.R