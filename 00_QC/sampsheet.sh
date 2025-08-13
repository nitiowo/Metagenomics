#!/bin/bash

cd ../raw_data
ls -1 *R1*.fastq | sed -E 's/_R1.*$//' | sort -u > samples.txt
