# alpha.R
# Alpha diversity: compute and plot metrics

library(phyloseq)
library(tidyverse)
library(ggplot2)

source("setup.R")
source("zoop_functions.R")

outdir <- file.path(output_root, "alpha")
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

# Control
use_ps_list <- ps_all_methods
use_markers <- NULL   # NULL = all markers
use_lakes   <- NULL   # NULL = all lakes
use_metrics <- c("Observed", "InvSimpson")
color_var   <- "Lake"
fill_colors <- lake_colors

# Filter
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)

# Compute a div across all markers
alpha_all <- compute_alpha_all(ps_filt, use_metrics, lake_order)

alpha_long <- alpha_all %>%
  pivot_longer(cols = all_of(use_metrics), names_to = "Metric", values_to = "Value")

# Plot a div faceted by marker
p_lk <- ggplot(alpha_long, aes(x = .data[[color_var]], y = Value,
                               fill = .data[[color_var]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_grid(Metric ~ Marker, scales = "free_y") +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Alpha diversity - between lakes (per marker)",
       y = "Diversity value") +
  theme(legend.position = "bottom")
save_plot(p_lk, file.path(outdir, "figures", "alpha_boxplot_lake_by_marker.pdf"),
          width = 16, height = 10)
