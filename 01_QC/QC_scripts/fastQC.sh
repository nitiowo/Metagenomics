#!/bin/bash --login

# Load fastqc (eg; HPC module or conda) if necessary
# conda activate fastqc-env

# Create output directory if it doesn't exist
mkdir -p ../fastqc_output

# Finding raw data
RAW_DIR=$(cat ../../00_prep/raw_dir.txt )

# Run fastqc on every file in raw data directory
fastqc "$RAW_DIR"/* -o ../fastqc_output/

# Deactivate fastqc
# conda deactivate if necessary
