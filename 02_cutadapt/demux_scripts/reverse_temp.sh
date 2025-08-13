#!/usr/bin/env bash --login
set -euo pipefail

# Usage: reverse_demux.sh <SAMPLE>
SAMPLE=$1
PRIMER_REV=../primers/reverse_primers.fa
IN_R1=../../raw_data/${SAMPLE}_R1_001.fastq.gz
IN_R2=../../raw_data/${SAMPLE}_R2_001.fastq.gz
OUTDIR=../demux_out/01a_reverse/
HCODIR=../demux_out/01a_reverse/HCO2198
SSUDIR=/demux_out/01a_reverse/SSUR22
UNKDIR=/demux_out/01a_reverse/unknown
LOGDIR=../demux_out/01a_reverse/logs
SUMDIR=../demux_out/summaries

mkdir -p "$OUTDIR" "$HCODIR" "$SSUDIR" "$LOGDIR" "$SUMDIR" 

# Helper functions to count reads and calculate average read length
count_reads(){ echo $(( $(zcat < "$1" | wc -l) / 4 )); }
avg_len() {
    zcat < "$1" | awk 'NR % 4 == 2 { total += length($0); count++ }
                     END { print total / count }'
}

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
  --info-file="${OUTDIR}/rev_info.tsv" \
  -p "${OUTDIR}/${SAMPLE}_{name}_R1.fastq.gz" \
  -o "${OUTDIR}/${SAMPLE}_{name}_R2.fastq.gz" \
  "${IN_R2}" "${IN_R1}" > $LOGDIR/${SAMPLE}_reverse.log 2>&1

### IMPORTANT NOTE!!!!!!
# During demultiplexing, cutadapt only uses forward primers to decide where each read is written.
# To 'cheat' this, we switch the positions of the reverse primers, reverse read input, and reverse read output
# So we 1) change -G to -g, 
# 2) switch -o and -p
# 3) switch the orders of the input files so the reverse reads come first.
# We do NOT do this in steps that use forward primers.

mv ${OUTDIR}/*HCO* $HCODIR
mv ${OUTDIR}/*SSU* $SSUDIR
mv ${OUTDIR}/*unk* $UNKDIR


# Count reads 
count_rev_hco_r=$(count_reads "${HCODIR}/${SAMPLE}_HCO2198_R2.fastq.gz" || echo 0)
count_rev_ssu_r=$(count_reads "${SSUDIR}/${SAMPLE}_SSUR22_R2.fastq.gz" || echo 0)
count_rev_unk_r=$(count_reads "${UNKDIR}/${SAMPLE}_unknown_R2.fastq.gz" || echo 0)

count_rev_hco_f=$(count_reads "${HCODIR}/${SAMPLE}_HCO2198_R1.fastq.gz" || echo 0)
count_rev_ssu_f=$(count_reads "${SSUDIR}/${SAMPLE}_SSUR22_R1.fastq.gz" || echo 0)
count_rev_unk_f=$(count_reads "${UNKDIR}/${SAMPLE}_unknown_R1.fastq.gz" || echo 0)

count_total_f=$(count_reads "")

# Calculate average read length
avglen_rev_hco_r=$(avg_len "${HCODIR}/${SAMPLE}_HCO2198_R2.fastq.gz" || echo 0)
avglen_rev_ssu_r=$(avg_len "${SSUDIR}/${SAMPLE}_SSUR22_R2.fastq.gz" || echo 0)
avglen_rev_unk_r=$(avg_len "${UNKDIR}/${SAMPLE}_unknown_R2.fastq.gz" || echo 0)

avglen_rev_hco_f=$(avg_len "${HCODIR}/${SAMPLE}_HCO2198_R1.fastq.gz" || echo 0)
avglen_rev_ssu_f=$(avg_len "${SSUDIR}/${SAMPLE}_SSUR22_R1.fastq.gz" || echo 0)
avglen_rev_unk_f=$(avg_len "${UNKDIR}/${SAMPLE}_unknown_R1.fastq.gz" || echo 0)

# If forward and reverse counts match, value is YES
if [[ $count_rev_hco_r == $count_rev_hco_f && \
  $count_rev_ssu_r == $count_rev_ssu_f && \
  $count_rev_unk_r == $count_rev_unk_f ]]; then
  countsMatch="YES"
else
  countsMatch="NO"
fi

# If the rev_count_summary file doesn't exist, create it and add column names.
if [ ! -f ${SUMDIR}/"rev_count_summary.tsv" ]; then
  echo -e "Sample\tCount_HCO\tCount_SSU\tCount_unk\tCountsMatch?" >> ${SUMDIR}/rev_count_summary.tsv
fi

# Add summary count into the rev_count_summary file
echo -e "${SAMPLE}\t${count_rev_hco_r}\t${count_rev_ssu_r}\t${count_rev_unk_r}\t${countsMatch}" >> ${SUMDIR}/rev_count_summary.tsv


# If the rev_readLen_summary file doesn't exist, create it and add column names.
if [ ! -f ${SUMDIR}/"rev_readLen_summary.tsv" ]; then
  echo -e "Sample\treadLen_HCO_R\treadLen_HCO_F \
  \treadLen_SSU_R\treadLen_SSU_F \
  \treadLen_unk_R\treadLen_unk_F" >> ${SUMDIR}/rev_readLen_summary.tsv
fi

# Add summary read lengths into the rev_readLen_summary file
echo -e "${SAMPLE}\t${avglen_rev_hco_r}\t${avglen_rev_hco_f} \
\t${avglen_rev_ssu_r}\t${avglen_rev_ssu_f} \
\t${avglen_rev_unk_r}\t${avglen_rev_unk_f}" >> ${SUMDIR}/rev_readLen_summary.tsv