# setup.R
# Load data, build phyloseq objects, define palettes, create combined datasets

library(phyloseq)
library(tidyverse)
library(vegan)
library(patchwork)
library(RColorBrewer)
library(flextable)
library(pheatmap)
library(viridis)
library(VennDiagram)

# Load config
if (!exists("pipeline_dir")) source("config.R")
setwd(pipeline_dir)
source("scripts/zoop_functions.R")

set.seed(42)

# ---- Palettes and constants ----
lake_order <- c("Superior", "Michigan", "Huron", "Erie", "Ontario")

lake_colors <- c(
  Superior = "#1b9e77", Michigan = "#d95f02", Huron = "#7570b3",
  Erie     = "#e7298a", Ontario  = "#66a61e"
)

marker_colors <- c(
  Folmer = "#e41a1c", Leray = "#377eb8",
  `18S`  = "#4daf4a", Morphology = "#984ea3"
)

# Standard plot settings
plot_base_size <- 12
plot_dpi       <- 300

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order",
               "Family", "Genus", "Species")

alpha_metrics <- c("Observed", "InvSimpson")

# ---- Load phyloseq objects ----
ps_folmer <- readRDS(rds_folmer)
ps_leray  <- readRDS(rds_leray)
ps_18S    <- readRDS(rds_18S)
ps_morph  <- readRDS(rds_morph)

# ---- Set lake ordering ----
ps_folmer <- set_lake_order(ps_folmer, lake_order)
ps_leray  <- set_lake_order(ps_leray, lake_order)
ps_18S    <- set_lake_order(ps_18S, lake_order)
ps_morph  <- set_lake_order(ps_morph, lake_order)

# ---- Named lists ----
ps_markers     <- list(Folmer = ps_folmer, Leray = ps_leray, `18S` = ps_18S)
ps_all_methods <- c(ps_markers, list(Morphology = ps_morph))
ps_coi         <- list(Folmer = ps_folmer, Leray = ps_leray)

# ---- P/A Versions ----
ps_markers_pa     <- lapply(ps_markers, to_pa)
ps_all_methods_pa <- lapply(ps_all_methods, to_pa)
ps_coi_pa         <- lapply(ps_coi, to_pa)

# ---- Station-Aggregated Objects ----
ps_morph_by_station <- ps_morph %>%
  rename_samples_to_station() %>%
  aggregate_to_station()

ps_markers_by_station <- lapply(ps_markers, function(ps) {
  ps %>% rename_samples_to_station() %>% aggregate_to_station()
})

# ---- Combined Datasets ----
shared_stations <- Reduce(intersect,
                          c(lapply(ps_markers_by_station, sample_names),
                            list(sample_names(ps_morph_by_station))))

ps_markers_combined <- combine_ps_pa(ps_markers_by_station, "Species",
                                      shared_stations)
ps_coi_combined <- combine_ps_pa(
  ps_markers_by_station[c("Folmer", "Leray")], "Species", shared_stations)
ps_all_combined <- combine_ps_pa(
  c(ps_markers_by_station, list(Morphology = ps_morph_by_station)),
  "Species", shared_stations)

# ---- Morph Info ----
morph_sample_names <- sample_names(ps_morph)
morph_station_ids  <- data.frame(sample_data(ps_morph))$Station_ID

# ---- Output Directories ----
analysis_dirs <- c("exploratory", "alpha", "beta", "composition",
                   "overlap", "heatmaps", "differential",
                   "geographic", "focal_taxon", "varpart", "trebitz_compare")
for (d in analysis_dirs) {
  dir.create(file.path(output_root, d, "figures"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_root, d, "stats"),
             recursive = TRUE, showWarnings = FALSE)
}
