#!/usr/bin/env bash
set -euo pipefail

# Usage: reverse_demux.sh <SAMPLE>
SAMPLE=$1
PRIMER_REV=primers/reverse_primers.fa
IN_R1=data/${SAMPLE}_R1.fastq.gz
IN_R2=data/${SAMPLE}_R2.fastq.gz
OUTDIR=demux/${SAMPLE}/reverse
mkdir -p "$OUTDIR"

# helper to count reads
count_reads(){ echo $(( $(zcat "$1" | wc -l) / 4 )); }

echo "Processing reverse demux for \$SAMPLE"

# run cutadapt, bin by reverse primers
cutadapt \
  -G file:"${PRIMER_REV}" \
  --pair-filter=any \
  --info-file="${OUTDIR}/rev_info.tsv" \
  -o "${OUTDIR}/${SAMPLE}_{name}_R1.fastq.gz" \
  -p "${OUTDIR}/${SAMPLE}_{name}_R2.fastq.gz" \
  "${IN_R1}" "${IN_R2}" > logs/${SAMPLE}_reverse.log 2>&1

# count reads
rev_hco=$(count_reads "${OUTDIR}/HCO2198_R1.fastq.gz"  || echo 0)
rev_ssu=$(count_reads "${OUTDIR}/SSUR_R1.fastq.gz"     || echo 0)
rev_unm=$(count_reads "${OUTDIR}/unmatched_R1.fastq.gz" || echo 0)

echo -e "${SAMPLE}\t${rev_hco}\t${rev_ssu}\t${rev_unm}" >> summaries/rev_summary.tsv