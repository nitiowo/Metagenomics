#!/bin/bash

# This script generates a sample sheet from your raw data directory

cd ../raw_data
ls -1 *R1*.fastq.gz | sed -E 's/_R1.*$//' | sort -u > ../00_prep/samples.txt
