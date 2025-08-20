# This script creates a bash script that runs a bunch of cat commands to concatenate reads together in specific ways.
# Usage: Rscript make_cat_script.R <total_summary_csv> <total_manifest_csv>

library(dplyr)

setwd('./')

args <- commandArgs(trailingOnly = TRUE)

# sum_df <- read.csv("../demux_out/01b_forward/summaries/total_fwd_out_summary.csv", header = T)
# man_df <- read.csv("../demux_out/01b_forward/summaries/total_fwd_out_manifest.csv", header = T)

sum_df  <- read.csv(args[1], header = T)
man_df <- read.csv(args[2], header = T)

prim_combos <- tibble(
  reverse_bin = c("HCO2198", "HCO2198", "SSUR22", "SSUR22", "unknown", "unknown", "unknown"),
  fwd_bin     = c("LCO1490", "mICOintF", "SSUF04", "unknown", "LCO1490", "mICOintF", "SSUF04"),
  pair_id =c(1,2,3,4,5,6,7)
)

paths_df <- sum_df %>% 
  inner_join(prim_combos, by = c("reverse_bin", "fwd_bin")) %>% 
  inner_join(man_df, by = c("sample","reverse_bin","fwd_bin")) %>%
  select(sample, reverse_bin, fwd_bin, out_R1_path, out_R2_path, pair_id)

cat_groups <- list(LCO1490=c(1,5), mICOintF=c(2,6), SSUF04=c(3,4,7))

samples <- unique(paths_df$sample)

results_list <- list()

for (s in samples){
  for (g in names(cat_groups)){
    curr_pairs <- cat_groups[[g]]
    
    samp_group <- paths_df %>%
      filter(pair_id %in% curr_pairs, sample == s)
    
    F_paths <- samp_group %>% pull(out_R1_path) %>% paste(collapse = " ")
    R_paths <- samp_group %>% pull(out_R2_path) %>% paste(collapse = " ")
    
    f_command <- paste0("cat ", F_paths, " > $1/", s, "_", g, "_", "R1.fastq.gz")
    r_command <- paste0("cat ", R_paths, " > $1/", s, "_", g, "_", "R2.fastq.gz")
    
    row1 <- data.frame(sample = s, group = g, f_r = "F", cat_command = f_command)
    row2 <- data.frame(sample = s, group = g, f_r = "R", cat_command = r_command)
    
    results_list[[length(results_list) + 1]] <- row1
    results_list[[length(results_list) + 1]] <- row2
  }
}

path_combos <- bind_rows(results_list)

writeLines(c("#!/bin/bash", path_combos$cat_command), "rescue_add.sh")
