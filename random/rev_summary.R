#!/usr/bin/env Rscript
# rev_summary.R
#
# Usage:
#   Rscript rev_summary.R manifest.csv outdir
#
# manifest.csv must contain columns:
#   sample, reverse_bin, out_R1_path, out_R2_path
#
# Output files (in outdir):
#   <SAMPLE>_rev_countLen.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(R.utils)
  library(Biostrings)
})

# Take 2 arguments as input
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript rev_summary.R manifest.csv outdir")
}

# Save the two arguments (manifest path and output directory) to variables
manifest_path <- args[1]
outdir <- args[2]

raw_path <- "../../raw_data/"

# Provide the following two lines manually when running in interactive mode or RStudio
# manifest_path <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/02_cutadapt/demux_out/01a_reverse/summaries/manifests/01_S01_L001_rev_out_man.csv"
# outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/02_cutadapt/demux_out/01a_reverse/summaries/summary_data/"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)  # Double-check that outdir exists

manifest <- read_csv(manifest_path, col_types = cols(
  CA_exit_status = col_logical()
))

# Creating empty results df
results <- manifest %>%
  select(-(out_R1_path:CA_exit_status)) %>% 
  mutate(
    Count_R = 0L,
    Count_F = 0L,
    AvgLen_R = 0.0,
    AvgLen_F = 0.0,
    FR_countsMatch = F,
    binSum_countsMatch = F
  )

message("Computing counts and average lengths for manifest rows")

# Issue warning if more than 1 sample name detected
if (length(unique(manifest$sample)) == 1) {
  samp <- manifest$sample[1]
} else {
  message("WARNING!!! More than 1 sample name detected!")
}

# Initialize count sums
countSum_f <- 0
countSum_r <- 0

# For each row (primer), calculate read counts and average read lengths, save them to results df
for (i in 1:nrow(manifest)) {
  r1 <- manifest$out_R1_path[i]   # Forward read path in your manifest (R1)
  r2 <- manifest$out_R2_path[i]   # Reverse read path (R2)
  r1_fq <- readDNAStringSet(r1, format = "fastq")  # Read in the whole fastqs
  r2_fq <- readDNAStringSet(r2, format = "fastq")
  results$Count_F[i] <- length(r1_fq)   # Calculate and save read counts and average lengths
  results$Count_R[i] <- length(r2_fq)
  results$AvgLen_F[i] <- mean(width(r1_fq))
  results$AvgLen_R[i] <- mean(width(r2_fq))
  
  if (results$Count_F[i] == results$Count_R[i]) {  # Check if F and R counts match for each primer
    results$FR_countsMatch <- TRUE
  } else {
    results$FR_countsMatch <- FALSE
  }
  if (results$reverse_bin[i] != "raw") {  # Update count totals (expect for the raw data read counts)
    countSum_f <- countSum_f + cf
    countSum_r <- countSum_r + cr
  }
}

rawcount_F <- results$Count_F[results$reverse_bin=="raw"]
rawcount_R <- results$Count_R[results$reverse_bin=="raw"]
if (countSum_f == rawcount_F & countSum_r == rawcount_R) {   # Check if primer count sums equal total raw data counts
  results$binSum_countsMatch[] <- TRUE
}
  
write_csv(results, file = paste0(outdir, samp, "_rev_countLen.csv"))