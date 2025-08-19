#!/bin/bash

# This script generates a sample sheet from your raw data directory

# Usage: bash sampsheet.sh

RAW_DIR=$(cat raw_dir.txt)

ls -1 "$RAW_DIR"/*R1*.fastq.gz \
  | xargs -n1 basename \
  | sed -E 's/_R1.*$//' \
  | sort -u \
 > ./samples.txt
