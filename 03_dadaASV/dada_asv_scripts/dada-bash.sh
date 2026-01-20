#!/bin/bash
#$ -M nvincen2@nd.edu
#$ -m abe
#$ -q largemem
#$ -pe smp 8
#$ -N dada2_mICOint
#$ -o job_logs/
#$ -e job_logs/

conda activate dada2-env

Rscript /temp180/mpfrende/nvincen2/Metagenomics/03_dadaASV/dada_asv_scripts/dada2_asv.R\
	--input /temp180/mpfrende/nvincen2/Metagenomics/02_cutadapt/demux_out/02c_rescue/mICOintF \
	--output /temp180/mpfrende/nvincen2/Metagenomics/03_dadaASV/dada_asv_output/mICOintF \
	--outprefix "mICOintF"
