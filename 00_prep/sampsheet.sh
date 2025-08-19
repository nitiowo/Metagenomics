#!/bin/bash

# This script generates a sample sheet from your raw data directory

# Usage: bash sampsheet.sh

RAW_DIR=$(cat raw_dir.txt)

cd $RAW_DIR
ls -1 *R1*.fastq.gz | sed -E 's/_R1.*$//' | sort -u > ../00_prep/samples.txt
