#!/bin/bash --login

# Load multiqc (through module or conda)
# conda activate fastqc-env

# Create output directory if it doesn't exist
mkdir -p ../multiqc_output

# Run multiqc on fastqc output
multiqc ../fastqc_output/ -o ../multiqc_output

# Deactivate multiqc
# conda deactivate
