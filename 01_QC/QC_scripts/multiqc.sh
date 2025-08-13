#!/bin/bash --login

# Load multiqc (eg; HPC module or conda) if necessary
#conda activate fastqc-env

# Create output directory if it doesn't exist
mkdir -p ../multiqc_out

# Run multiqc on fastqc output
multiqc ../fastqc_out/ -o ../multiqc_out

# Deactivate multiqc
#conda deactivate if necessary
