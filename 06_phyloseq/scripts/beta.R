# beta.R
# Beta diversity: ordination

library(phyloseq)
library(tidyverse)
library(vegan)
library(ggplot2)

source("setup.R")

outdir <- file.path(output_root, "beta")
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

# Control
use_ps      <- ps_all_methods[["COI"]]  
color_var   <- "Lake"
fill_colors <- lake_colors
ord_method  <- "NMDS"
dist_method <- "jaccard"
use_binary  <- TRUE
plot_title  <- "NMDS Jaccard (P/A) - COI"
out_file    <- file.path(outdir, "figures", "nmds_jaccard.pdf")
plot_width  <- 8
plot_height <- 6

# Presence/absence transform
if (use_binary) {
  otu <- as(otu_table(use_ps), "matrix")
  otu[otu > 0] <- 1
  otu_table(use_ps) <- otu_table(otu, taxa_are_rows = taxa_are_rows(use_ps))
}

# Ordination
ord <- ordinate(use_ps, method = ord_method, distance = dist_method)

# Plot
p <- plot_ordination(use_ps, ord, color = color_var) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = fill_colors) +
  theme_minimal(base_size = 10) +
  labs(title = plot_title)

ggsave(out_file, p, width = plot_width, height = plot_height)
