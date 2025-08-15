#!/usr/bin/env bash --login

# This script calculates counts and read length averages for output files from cutadapt steps.
# This process has since been moved to an R script, so this script is deprecated.


set -euo pipefail # Kills script in case of command failure


# _________________ Summarizing outputs _____________________


primers=("HCO2198" "SSUR22" "unknown")      # Order matters for summaries
dirs=("$HCODIR" "$SSUDIR" "$UNKDIR")
short=("HCO" "SSU" "unk")                   # Short names used in headers/output

# Create associative arrays to store values (using strings as keys)
declare -A count_rev_r count_rev_f avglen_rev_r avglen_rev_f

# Initalize counts. Set countsMatch to YES until we see a mismatch
count_rev_total_r=0
count_rev_total_f=0
FR_countsMatch="YES"

# Single loop: compute counts, avg lengths, totals, and FR match per each primer
for i in "${!primers[@]}"; do  # ${!primers[@]} expands to indeces of all elements in array
  primer="${primers[i]}" # Save primer in first index position
  dir="${dirs[i]}"

  rpath="${dir}/${SAMPLE}_${primer}_R2.fastq.gz"
  fpath="${dir}/${SAMPLE}_${primer}_R1.fastq.gz"

  # Compute counts and averages
  count_r=$(count_reads "$rpath" 2>/dev/null || echo 0)
  count_f=$(count_reads "$fpath" 2>/dev/null || echo 0)
  avg_r=$(avg_len "$rpath" 2>/dev/null || echo 0)
  avg_f=$(avg_len "$fpath" 2>/dev/null || echo 0)

  # Save current count to appropriately appropriate array and primer index
  count_rev_r["$primer"]=$count_r
  count_rev_f["$primer"]=$count_f
  avglen_rev_r["$primer"]=$avg_r
  avglen_rev_f["$primer"]=$

  # Total counts
  count_rev_total_r=$((count_rev_total_r + count_r))
  count_rev_total_f=$((count_rev_total_f + count_f))

  # Check forward/reverse per-primer; if any mismatch, flag as NO
  if [[ $count_r -ne $count_f ]]; then
    FR_countsMatch="NO"
  fi
done

# Raw data read counts
count_raw_r=$(count_reads "$IN_R2" 2>/dev/null || echo 0)
count_raw_f=$(count_reads "$IN_R1" 2>/dev/null || echo 0)

# Check if read totals match
if [[ $count_rev_total_r -eq $count_raw_r && $count_rev_total_f -eq $count_raw_f ]]; then
  Tot_countsMatch="YES"
else
  Tot_countsMatch="NO"
fi

# Helper to ensure header exists
# Takes two arguments: file to write to, and header info
ensure_file_with_header() {
  local file="$1"; local header="$2"
  if [[ ! -f "$file" ]]; then
    printf "%s\n" "$header" > "$file"
  fi
}

# Write count summary file
count_file="${SUMDIR}/rev_count_summary.tsv"
count_header="Sample"                   # Create first column
for s in "${short[@]}"; do            # For each item in "short" primer names list:
    count_header+=$'\t'"Count_${s}"     # Create new column
done
count_header+=$'\t'"FR_CountsMatch?"$'\t'"Tot_CountsMatch?"  # Add countsMatch columns
ensure_file_with_header "$count_file" "$count_header"        # Appends header line to summary file

# Build each line by adding values to each column in primer order
line="$SAMPLE"                            # Value in first column
for primer in "${primers[@]}"; do       # For each primer:
    line+=$'\t'"${count_rev_r[$primer]}"  # Add count to the line (tab separated)
done
line+=$'\t'"${FR_countsMatch}"$'\t'"${Tot_countsMatch}"   # Add countsMatch values

printf "%s\n" "$line" >> "$count_file"

# Write read length summary file
len_file="${SUMDIR}/rev_readLen_summary.tsv"
len_header="Sample"
for s in "${short[@]}"; do
    len_header+=$'\t'"readLen_${s}_R"$'\t'"readLen_${s}_F"
done
ensure_file_with_header "$len_file" "$len_header"

# Build lines for read lengths
len_line="$SAMPLE"
for primer in "${primers[@]}"; do
    len_line+=$'\t'"${avglen_rev_r[$primer]}"$'\t'"${avglen_rev_f[$primer]}"
done

printf "%s\n" "$len_line" >> "$len_file"