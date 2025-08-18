#!/usr/bin/env bash --login

# Usage: bash total_summaries.sh

PRIMER_REV=../primers/reverse_primers.fa

TOTAL_MAN=../demux_out/01a_reverse/summaries/total_manifest.csv
TOTAL_SUM=../demux_out/01a_reverse/summaries/total_summary.csv

# Extract primer names from fasta file
mapfile -t PRIMERS < <(grep '^>' "$PRIMER_REV" | sed 's/^>//' | tr -d '\r')
# Append unknown
PRIMERS+=("unknown")

num_prim=${#PRIMERS[@]}    # Number of primers
num_lines=$((num_prim + 1))  # Number of lines to extract (no. of primers + raw)

man_dir=../demux_out/01a_reverse/summaries/manifests/
man_files=../demux_out/01a_reverse/summaries/manifests/*
first_man=$(ls "$man_dir" | head -n 1)    # Grabs first file in directory
head -n 1 "$man_dir/$first_man" > "$TOTAL_MAN"   # Grabs first line from that file (header line)

echo "Creating total manifest"
for file in ${man_files[@]}; do
	tail -n $num_lines $file >> $TOTAL_MAN  # Grabs last n lines from manifest file (n = number of bins)
done

sum_dir=../demux_out/01a_reverse/summaries/summary_data/
sum_files=../demux_out/01a_reverse/summaries/summary_data/*
first_sum=$(ls "$sum_dir" | head -n 1)
head -n 1 "$sum_dir/$first_sum" > "$TOTAL_SUM"

echo "Creating total summary"
for file in ${sum_files[@]}; do
	tail -n $num_lines $file >> $TOTAL_SUM
done