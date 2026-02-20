# setup.R
# Load all data, build phyloseq objects, define palettes

library(phyloseq)

#setwd('/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/')
set.seed(42)
output_root <- "output"

# Palettes and constants --    
lake_order <- c("Superior", "Michigan", "Huron", "Erie", "Ontario")

lake_colors <- c(
  Superior = "#1b9e77", Michigan = "#d95f02", Huron = "#7570b3",
  Erie     = "#e7298a", Ontario  = "#66a61e"
)

marker_colors <- c(
  Folmer = "#e41a1c", Leray = "#377eb8",
  `18S`  = "#4daf4a", Morphology = "#984ea3"
)

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order",
               "Family", "Genus", "Species")

alpha_metrics <- c("Observed", "InvSimpson")

# Load phyloseq objects 
ps_folmer <- readRDS("setup_output/folmer_ps.RDS")
ps_leray  <- readRDS("setup_output/leray_ps.RDS")
ps_18S    <- readRDS("setup_output/ssu_ps.RDS")
ps_morph  <- readRDS("setup_output/morph_ps.RDS")

#  Named lists 
ps_markers     <- list(Folmer = ps_folmer, Leray = ps_leray, `18S` = ps_18S)
ps_all_methods <- c(ps_markers, list(Morphology = ps_morph))
ps_coi         <- list(Folmer = ps_folmer, Leray = ps_leray)

#  Morph info 
morph_sample_names <- sample_names(ps_morph)
morph_station_ids  <- data.frame(sample_data(ps_morph))$Station_ID