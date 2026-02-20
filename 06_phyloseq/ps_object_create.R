# This script loads dada2 output, standardizes metadata and phyloseq formats, and outputs phyloseq objects

library(phyloseq)
library(dplyr)
library(Biostrings)

setwd("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/")
dir.create("06_phyloseq/setup_output")
outdir <- ("06_phyloseq/setup_output/")

folmer_taxa <- readRDS("04_dadaTax/dada_tax_output/LCO1490/Folmer_taxa.Rds")
folmer_seqtab <- readRDS("03_dadaASV/dada_asv_output/LCO1490/LCO1490_seqtab.nochim.Rds")
folmer_taxtable <- tax_table(folmer_taxa)
folmer_otu <- otu_table(folmer_seqtab, taxa_are_rows = FALSE)

leray_taxa <- readRDS("04_dadaTax/dada_tax_output/mICOintF/Leray_taxa.Rds")
leray_seqtab <- readRDS("03_dadaASV/dada_asv_output/mICOintF/mICOintF_seqtab.nochim.Rds")
leray_taxtable <- tax_table(leray_taxa)
leray_otu <- otu_table(leray_seqtab, taxa_are_rows = FALSE)

ssu_taxa <- readRDS("04_dadaTax/dada_tax_output/SSUF04/SSU_taxa.Rds")
ssu_seqtab <- readRDS("03_dadaASV/dada_asv_output/SSUF04/SSU_seqtab.nochim.Rds")
ssu_taxtable <- tax_table(ssu_taxa)
ssu_otu <- otu_table(ssu_seqtab, taxa_are_rows = FALSE)

morph_taxa <- read.csv("data/morph_data/morph_ps_taxa.csv", header = TRUE, row.names = 1)
morph_taxtable <- tax_table(as.matrix(morph_taxa))

morph_counts <- read.csv("data/morph_data/morph_ps_otu.csv", header = TRUE, check.names = FALSE, 
                         na.strings = c("NA", ""), row.names = 1)

# Replace NAs with 0s and convert all numeric columns to integers
morph_counts <- morph_counts %>%
  mutate(across(everything(), ~ replace(., is.na(.), 0))) %>%
  mutate(across(where(is.numeric), ~ as.integer(round(.))))

morph_otu <- otu_table(morph_counts, taxa_are_rows = FALSE)

metadata <- read.csv("data/zoop96_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
metadata$Mesh <- as.factor(metadata$Mesh)
meta_samples <- sample_data(metadata)
sample_names(meta_samples) <- metadata$Sample_ID

morph_meta <- read.csv("data/morph_data/morph_ps_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
morph_samples <- sample_data(morph_meta)
sample_names(morph_samples) <- morph_meta$Sample_ID

ps_folmer <- phyloseq(folmer_otu, folmer_taxtable, meta_samples)
ps_leray <- phyloseq(leray_otu, leray_taxtable, meta_samples)
ps_ssu <- phyloseq(ssu_otu, ssu_taxtable, meta_samples)
ps_morph <- phyloseq(morph_otu, morph_taxtable, morph_samples)

# Rename and save
ps_list <- list(folmer = ps_folmer, leray = ps_leray, ssu = ps_ssu, morph = ps_morph)

for (nm in names(ps_list)) {
  fileprefix <- paste0(outdir, nm)
  
  res <- rename_asvs_from_empty_refseq(
    ps_list[[nm]],
    save_map   = paste0(fileprefix, "_asv_map.csv"),
    save_fasta = paste0(fileprefix, "_asv.fasta")
  )
  
  # Replace object inside the list, save RDS objects
  ps_list[[nm]] <- res$ps
  saveRDS(ps_list[[nm]], file = paste0(fileprefix, "_ps.RDS"))
}