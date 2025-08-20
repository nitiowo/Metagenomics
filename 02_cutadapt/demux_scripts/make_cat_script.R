# This script creates a bash script that runs a bunch of cat commands to concatenate reads together in specific ways.
# Usage: Rscript make_cat_script.R <total_manifest_csv>

library(dplyr)

setwd('./')

args <- commandArgs(trailingOnly = TRUE)

# Stop script if not enough args supplied
if (length(args) < 1) {
  stop("Usage: Rscript make_cat_script.R <total_manifest_csv>")
}

# Uncomment the following  line to run interactively
# man_df <- read.csv("../demux_out/01b_forward/summaries/total_fwd_out_manifest.csv", header = T)

man_df <- read.csv(args[1], header = T)

# Checking for required columns
stopifnot(all(c("sample","reverse_bin","fwd_bin") %in% names(man_df)))

# Setting primer combos of interest
# These sets are not being added together - these are unique combinations of reverse bin + forward bin
prim_combos <- tibble(
  reverse_bin = c("HCO2198", "HCO2198", "SSUR22", "SSUR22", "unknown", "unknown", "unknown"),  # List names must match column names in mandf and sumdf
  fwd_bin     = c("LCO1490", "mICOintF", "SSUF04", "unknown", "LCO1490", "mICOintF", "SSUF04"),
  pair_id =c(1,2,3,4,5,6,7)   # Giving each pair a unique ID
)

# Subsetting big dataframe
paths_df <- man_df %>% 
  # Joining man_df with prim_combos (only selected combos will remain) - add pair_id column
  inner_join(prim_combos, by = c("reverse_bin", "fwd_bin")) %>% 
  select(sample, reverse_bin, fwd_bin, out_R1_path, out_R2_path, pair_id)

# Inspect paths_df and select how the pairs will be added together. Use the unique pair_ids.
# Give each group a name - these will go in output file names
cat_groups <- list(LCO1490=c(1,5), mICOintF=c(2,6), SSUF04=c(3,4,7))

samples <- unique(paths_df$sample)
results_list <- list()

for (s in samples){  # s takes on sample names
  for (g in names(cat_groups)){ # g takes on group names provided above
    curr_pairs <- cat_groups[[g]]   # Grab numbers inside the group
    
    samp_group <- paths_df %>%
      filter(pair_id %in% curr_pairs, sample == s)  # Grab only rows for that sample, and the pair IDs in the group
    
    F_paths <- samp_group %>% pull(out_R1_path) %>% paste(collapse = " ")  # Collapse all items the out_R1_path column, separated by spaces
    R_paths <- samp_group %>% pull(out_R2_path) %>% paste(collapse = " ")  # R2
    
    # Build cat command - uses $1/ so that you can provide output directory path as argument to the cat script
    f_command <- paste0("cat ", F_paths, " > $1/", s, "_", g, "_", "R1.fastq.gz")
    r_command <- paste0("cat ", R_paths, " > $1/", s, "_", g, "_", "R2.fastq.gz")
    
    row1 <- data.frame(sample = s, group = g, f_r = "F", cat_command = f_command)  # Build each row in output df
    row2 <- data.frame(sample = s, group = g, f_r = "R", cat_command = r_command)
    
    results_list[[length(results_list) + 1]] <- row1  # Add forward cat command to list
    results_list[[length(results_list) + 1]] <- row2  # Reverse
  }
}

# Combine all rows in results_list into a df
path_combos <- bind_rows(results_list)

# Prepend shebang and write each value in cat_command column in path_combos to new file
writeLines(c("#!/bin/bash", path_combos$cat_command), "rescue_add.sh")
