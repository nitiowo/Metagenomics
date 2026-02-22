# exploratory.R
# Dataset summaries

library(phyloseq)
library(tidyverse)

source("setup.R")
source("zoop_functions.R")

outdir <- file.path(output_root, "exploratory")
dir.create(file.path(outdir, "stats"), recursive = TRUE, showWarnings = FALSE)

# Control
use_ps      <- ps_all_methods[["Folmer"]]
use_name    <- "Folmer"
use_ranks   <- tax_ranks

# Summarise
summary_df <- summarise_ps(use_ps, use_name, use_ranks)

summary_df <- summary_df %>%
  select(dataset, samples, total_asvs, rank, unique_taxa,
         unassigned_asvs, pct_reads_unassigned)

write.csv(summary_df,
          file.path(outdir, "stats", "dataset_summary.csv"),
          row.names = FALSE)
