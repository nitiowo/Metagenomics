# alpha.R
# Alpha diversity: compute and plot metrics

library(phyloseq)
library(tidyverse)
library(ggplot2)

source("setup.R")

outdir <- file.path(output_root, "alpha")
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

# Control
use_ps      <- ps_all_methods[["COI"]]
color_var   <- "Lake"
fill_colors <- lake_colors
use_metrics <- c("Observed", "InvSimpson") # Shannon not needed
plot_title  <- "Alpha diversity - 18S conversion"
out_file    <- file.path(outdir, "figures", "alpha_boxplot.pdf")
plot_width  <- 10
plot_height <- 6

#  Compute alpha diversity
alpha_df <- estimate_richness(use_ps, measures = use_metrics)
alpha_df$Sample_ID <- sample_names(use_ps)

meta <- data.frame(sample_data(use_ps))
meta$Sample_ID <- rownames(meta)

alpha_df <- left_join(alpha_df, meta, by = "Sample_ID")

alpha_long <- alpha_df %>%
  pivot_longer(cols = all_of(use_metrics), names_to = "Metric", values_to = "Value")

# Plot
p <- ggplot(alpha_long, aes(x = .data[[color_var]], y = Value, fill = .data[[color_var]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 10) +
  labs(title = plot_title, y = "Diversity value") +
  theme(legend.position = "bottom")

ggsave(out_file, p, width = plot_width, height = plot_height)
