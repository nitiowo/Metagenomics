# This script loads dada2 output, standardizes metadata and phyloseq formats, and outputs phyloseq objects

library(phyloseq)
library(dplyr)

setwd("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/")
folmer_taxa <- readRDS("04_dadaTax/dada_tax_output/LCO1490/Folmer_taxa.Rds")
folmer_taxtable <- tax_table(folmer_taxa)

leray_taxa <- readRDS("04_dadaTax/dada_tax_output/mICOintF/Leray_taxa.Rds")
leray_taxtable <- tax_table(leray_taxa)

ssu_taxa <- readRDS("04_dadaTax/dada_tax_output/SSUF04/SSU_taxa.Rds")
ssu_taxtable <- tax_table(ssu_taxa)

folmer_seqtab <- readRDS("03_dadaASV/dada_asv_output/LCO1490/LCO1490_seqtab.nochim.Rds")
folmer_otu <- otu_table(folmer_seqtab, taxa_are_rows = FALSE)

leray_seqtab <- readRDS("03_dadaASV/dada_asv_output/mICOintF/mICOintF_seqtab.nochim.Rds")
leray_otu <- otu_table(leray_seqtab, taxa_are_rows = FALSE)

ssu_seqtab <- readRDS("03_dadaASV/dada_asv_output/SSUF04/SSU_seqtab.nochim.Rds")
ssu_otu <- otu_table(ssu_seqtab, taxa_are_rows = FALSE)

metadata <- read.csv("06_phyloseq/zoop96_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
metadata$Mesh <- as.factor(metadata$Mesh)
meta_samples <- sample_data(metadata)
sample_names(meta_samples) <- metadata$Sample_ID

ps_folmer <- phyloseq(folmer_otu, folmer_taxtable, meta_samples)
ps_leray <- phyloseq(leray_otu, leray_taxtable, meta_samples)
ps_ssu <- phyloseq(ssu_otu, ssu_taxtable, meta_samples)

morph_otu <- read.csv("Data/MORPH_counts_biomass.csv", header = TRUE, check.names = FALSE, na.strings = c("NA", ""))

# Replace NAs with 0s and convert all numeric columns to integers
morph_otu <- morph_otu %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(across(where(is.numeric), as.integer))



