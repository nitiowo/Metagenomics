#!/usr/bin/env bash
# Usage: samplesheet_gen.sh <inputDir> <outputFile>

# Makes sample sheet from fastq filenames

INDIR="${1:-.}"
OUT="${2:-samplesheet.tsv}"

printf "sampleID\tforwardReads\treverseReads\n" > "$OUT"

shopt -s nullglob
for f in "$INDIR"/*_R1*.fastq.gz; do
  [ -f "$f" ] || continue
  base=$(basename "$f")
  sample="${base%%_*}"
  sample="Z${sample}"  # Prefix Z to avoid numeric-only sample names

  # Find matching R2
  rev="${f/_R1/_R2}"
  if [ ! -f "$rev" ]; then
    rev="${f/R1/R2}"
  fi

  if [ -f "$rev" ]; then
    rev_abs=$(realpath "$rev")
  else
    rev_abs=""
  fi

  f_abs=$(realpath "$f")
  printf "%s\t%s\t%s\n" "$sample" "$f_abs" "$rev_abs" >> "$OUT"
done

echo "Wrote $OUT"
