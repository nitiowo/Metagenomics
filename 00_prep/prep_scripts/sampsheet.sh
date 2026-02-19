#!/bin/bash

# This script generates a sample sheet from your raw data directory
# Usage: bash sampsheet.sh

# Source config (create `config.sh` from `config.sh.example`)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

ls -1 "$RAW_DIR"/*R1*.fastq.gz \
  | xargs -n1 basename \
  | sed -E 's/_R1.*$//' \
  | sort -u \
 > "$SAMPLES_FILE"
