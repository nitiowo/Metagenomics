#!/usr/bin/env bash
# usage: samplesheet_gen.sh <inputDir> <outputFile>

INDIR="${1:-.}"
OUT="${2:-samplesheet.tsv}"

# small helper to get absolute path
abspath() {
  if command -v realpath >/dev/null 2>&1; then
    realpath -m "$1"
  elif command -v readlink >/dev/null 2>&1; then
    readlink -f "$1"
  else
    python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$1"
  fi
}

# header
printf "sampleID\tforwardReads\treverseReads\n" > "$OUT"

# iterate forward reads (matches examples like *_R1*.fastq.gz)
shopt -s nullglob
for f in "$INDIR"/*_R1*.fastq.gz; do
  [ -f "$f" ] || continue
  base=$(basename "$f")
  sample="${base%%_*}"                    # first token before first underscore
  sample="Z${sample}"		# add "Z" before sample names to avoid issues with numerical fields. Comment out if unneeded.
  # expected reverse file by simple substitution
  rev="${f/_R1/_R2}"
  # if that doesn't exist try also replacing "R1"->"R2" without underscore
  if [ ! -f "$rev" ]; then
    rev="${f/R1/R2}"
  fi
  # if still missing, leave empty
  if [ -f "$rev" ]; then
    rev_abs=$(abspath "$rev")
  else
    rev_abs=""
  fi
  f_abs=$(abspath "$f")
  printf "%s\t%s\t%s\n" "$sample" "$f_abs" "$rev_abs" >> "$OUT"
done

echo "Wrote $OUT"
