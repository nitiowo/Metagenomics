# composition.R
# Stacked barplots of relative abundance

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

source("setup.R")
source("zoop_functions.R")

outdir <- file.path(output_root, "composition")
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

# Control
use_ps     <- ps_all_methods[["Folmer"]]
use_rank   <- "Phylum"
top_n      <- 15
plot_title <- "Relative abundance - Folmer Phylum"
out_file   <- file.path(outdir, "figures", "barplot_phylum_folmer.pdf")
plot_width  <- 12
plot_height <- 8

# Aggregate to rank and get relative abundance
ps_agg <- agg_rank(use_ps, use_rank)
ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))

df <- psmelt(ps_rel)
df$Lake <- factor(df$Lake, levels = lake_order, ordered = TRUE)

# Top N taxa, lump everything else as "Other"
top_taxa <- df %>%
  group_by(.data[[use_rank]]) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = top_n) %>%
  pull(.data[[use_rank]])

df[[use_rank]] <- if_else(df[[use_rank]] %in% top_taxa,
                          as.character(df[[use_rank]]), "Other")
df[[use_rank]] <- factor(df[[use_rank]],
                         levels = c(as.character(top_taxa), "Other"))

# Make colors
fill_pal <- c(colorRampPalette(brewer.pal(12, "Set3"))(top_n), "grey70") # Set 3 has the most?
names(fill_pal) <- levels(df[[use_rank]])

# Plot
p <- ggplot(df, aes(x = Sample, y = Abundance, fill = .data[[use_rank]])) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = fill_pal) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        legend.position = "right") +
  labs(title = plot_title, y = "Relative abundance", x = "")

ggsave(out_file, p, width = plot_width, height = plot_height)
