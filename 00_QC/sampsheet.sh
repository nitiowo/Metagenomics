#!/bin/bash --login

ls -1 ../raw_data/*_R1_*.fastq | sed -E 's/_S[0-9]+_L[0-9]+_R1_001.fastq$//' | sort -u > ../raw_data/samples.txt

