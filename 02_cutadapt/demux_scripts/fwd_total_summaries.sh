#!/usr/bin/env bash --login

# Usage: bash total_summaries.sh

# Input directories
PRIMER_FWD=../primers/fwd_primers.fa
man_dir=../demux_out/01b_forward/summaries/manifests/
sum_dir=../demux_out/01b_forward/summaries/summary_data/

# Output files
TOTAL_MAN=../demux_out/01b_forward/summaries/total_fwd_out_manifest.csv
TOTAL_SUM=../demux_out/01b_forward/summaries/total_fwd_out_summary.csv

# Extract primer names from fasta file
mapfile -t PRIMERS < <(grep '^>' "$PRIMER_FWD" | sed 's/^>//' | tr -d '\r')
# Append unknown
PRIMERS+=("unknown")

# Getting primer info
num_prim=${#PRIMERS[@]}    # Number of primers
num_lines=$((num_prim + 1))  # Number of lines to extract (no. of primers + reverse parent) (raw line added with header)

# Comma-separated manifests; prints header once and only the first 'raw' line per sample
# Print header once
# Print raw only the first time per sample
# Otherwise print line
awk -F',' '
  FNR==1 && NR==1 { print; next }
  FNR==1 { next }                
  $2=="raw" { 
  	if (!seen[$1]++) print; next }      
  { print }                                       
' "$man_dir"/* > ${TOTAL_MAN}

# Creating manifest without raw info for forward array input
awk -F',' '$2 !~ /raw/ && $2 !~ /parent/' "${TOTAL_MAN}" > ../demux_out/01b_forward/summaries/total_fwd_noRawParent_manifest.csv

# Summary file
awk -F',' '
  FNR==1 && NR==1 { print; next }  
  FNR==1 { next }                
  $2=="raw" { if (!seen[$1]++) print; next }      
  { print }                                       
' "$sum_dir"/* > ${TOTAL_SUM}