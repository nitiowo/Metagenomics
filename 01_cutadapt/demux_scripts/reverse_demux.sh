#!/usr/bin/env bash
set -euo pipefail

# Usage: reverse_demux.sh <SAMPLE>

# Make sure cutadapt is in $PATH (module load or conda)

SAMPLE=$1
DATA_DIR=../raw_data
PRIMER_REV=./primers/reverse_primers.fa
IN_R1=${DATA_DIR}/${SAMPLE}_R1.fastq.gz
IN_R2=${DATA_DIR}/${SAMPLE}_R2.fastq.gz
OUTDIR=./cutadapt_out/${SAMPLE}/reverse

mkdir -p "$OUTDIR"
mkdir -p "${OUTDIR}/summaries"

# helper to count reads
ncount_reads(){ echo $(( $(zcat "$1" | wc -l) / 4 )); }

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
total_reads=$(count_reads "${DATA_DIR}/${SAMPLE}_R2.fastq.gz"  || echo 0)
rev_hco=$(count_reads "${OUTDIR}/${SAMPLE}_HCO2198_R1.fastq.gz"  || echo 0)
rev_ssu=$(count_reads "${OUTDIR}/${SAMPLE}_SSUR_R22.fastq.gz"     || echo 0)
rev_unm=$(count_reads "${OUTDIR}/${SAMPLE}_unmatched_R1.fastq.gz" || echo 0)

echo -e "Sample"\t"HCO_reads"\t"SSUR_reads"\t"unm_reads"\t"total_reads" > ${OUTDIR}/summaries/rev_summary.tsv
echo -e "${SAMPLE}\t${rev_hco}\t${rev_ssu}\t${rev_unm}\t$total_reads" >> ${OUTDIR}/summaries/rev_summary.tsv