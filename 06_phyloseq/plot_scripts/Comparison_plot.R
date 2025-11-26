library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library('ape')
library("vegan")
library("DESeq2")
library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
library(dendextend) ; packageVersion("dendextend") # 1.10.0
library(viridis) ; packageVersion("viridis") # 0.5.1
library(decontam)
library(dplyr)
library(Biostrings); packageVersion("Biostrings")
library("ggpubr")
library(scales)
library(microViz)
library(ggtext)
library(ggraph)
library(DT)
library(corncob)
library(shiny)
library(hash)

setwd('/Volumes/Samsung_1TB/Zooplankton/')
outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/taxonomic_distribution/"

####################################################################

# Importing data ----

####################################################################

SampleMetadata <- read.csv('/Volumes/Samsung_1TB/Zooplankton/Metagenomics/raw_data/zoop96_metadata.csv', header = TRUE, stringsAsFactors = TRUE)
seqtab.nochim <- readRDS("/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_SSUF/seqtab.nochim.Rds")
tax_info <- as.matrix(readRDS("/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_SSUF/taxa.Rds"))

marker <- "Leray"
markerdir <- paste0(outdir, marker)
dir.create(markerdir, recursive = TRUE)

# abundance <- readRDS('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_SSUF/abundance_SSUF.Rds')
# asv_ref <- read.csv('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_LCO/ASV_ref.csv')
# asv_fasta <- read.delim('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_LCO/ASVs.fa')

########### For Morphological ----

seqtab.nochim <- read.csv("/Volumes/Samsung_1TB/Zooplankton/phyloseq/morphological/morph_asv_samplesFixed.csv", header = TRUE)
tax_info <- as.matrix(readRDS("/Volumes/Samsung_1TB/Zooplankton/phyloseq/morphological/morph_taxClean.RDS"))

rownames(seqtab.nochim) <- seqtab.nochim$X
seqtab.nochim <- seqtab.nochim[,-1]

samp <- as.data.frame(colnames(seqtab.nochim))  # Creating metadata for morph - 2 columns - station_ID and Lake. 
colnames(samp) <- c("Station_ID") # (You could add year and month later)
samp <- samp %>%
  mutate(Lake = case_when(
    grepl("^HU", Station_ID) ~ "Huron",
    grepl("^ER", Station_ID) ~ "Erie",
    grepl("^SU", Station_ID) ~ "Superior",
    grepl("^MI", Station_ID) ~ "Michigan",
    grepl("^ON", Station_ID) ~ "Ontario",
    TRUE ~ NA_character_
  ))

rownames(samp) <- samp$Station_ID

seqtab.nochim[is.na(seqtab.nochim)] <- 0 # Replace NA with 0

seqtab.nochim <- seqtab.nochim %>%
  mutate(across(where(is.numeric), as.integer)) # Converts all numeric columns

OTU = otu_table(seqtab.nochim, taxa_are_rows = TRUE)
TAX = tax_table(tax_info)
SAM = sample_data(samp)

ps <- phyloseq(OTU, SAM, TAX)

saveRDS(ps, file = "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/morph_ps_integer.RDS")

# ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/phyloseq/morphological/morph_ps.RDS")

####################################################################

# Creating ps objects (skip to loading if already created) ----

####################################################################

rownames(SampleMetadata) <- rownames(seqtab.nochim)
SampleMetadata$Mesh <- as.factor(SampleMetadata$Mesh)
SampleMetadata$Year_Sampled <- as.factor(SampleMetadata$Year_Sampled)
SampleMetadata$Tube <- as.factor(SampleMetadata$Tube)

# Creating phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(SampleMetadata), 
               tax_table(tax_info))

sample_names(ps) <- paste0("S-", sample_names(ps))

# Rename ASV sequences to headers
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
ps

saveRDS(ps, file = paste0("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/", marker, "_dadaOriginal.RDS"))

###############################################################

############ For 18S ----
tax_info1 <- as.matrix(readRDS("/Volumes/Samsung_1TB/Zooplankton/phyloseq/SSUF/wrangling/SSU_animalia_taxFix_taxClean.RDS"))

ps_temp <- prune_taxa(rownames(tax_info1), ps)
ps_temp1 <- phyloseq(otu_table(otu_table(ps_temp)), 
                     tax_table(tax_info1),
                     refseq(refseq(ps_temp)),
                     sample_data(sample_data(ps_temp)))
ps <- ps_temp1
saveRDS(ps, file = paste0("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/", marker, "_dadaOriginal.RDS"))

####################################################################

# Load ps object if already created ----

####################################################################

Leray_ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/Leray_dadaOriginal.RDS")
Folmer_ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/Folmer_dadaOriginal.RDS")
ps_18S <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/18S_dadaOriginal.RDS")
morph_ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/morph_ps_integer.RDS")


outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/comparison_plots/"

# marker <- "Leray" # Only works with no slash at the end
# markerdir <- paste0(outdir, marker)
# dir.create(markerdir, recursive = TRUE)


random_tree <- rtree(ntaxa(orig_ps), rooted=TRUE, 
                     tip.label=taxa_names(orig_ps))
ps <- merge_phyloseq(orig_ps, random_tree) # Add tree to ps object



dim(otu_table(Leray_ps))
sample_sums(Leray_ps)

# Transform to incidence data (presence/absence)
Leray_inc_ps <- transform_sample_counts(Leray_ps, function(x) ifelse(x > 0, 1L, 0L))
Folmer_inc_ps <- transform_sample_counts(Folmer_ps, function(x) ifelse(x > 0, 1L, 0L))
ps_inc_18S <- transform_sample_counts(ps_18S, function(x) ifelse(x > 0, 1L, 0L))
morph_inc_ps <- transform_sample_counts(morph_ps, function(x) ifelse(x > 0, 1L, 0L))

sample_sums(tax_glom(Leray_inc_ps, taxrank = "Phylum"))

get_taxa_unique(Leray_ps, taxonomic.rank = "Class")

taxa_are_rows(Leray_inc_ps)

hello <- psmelt(Leray_ps)
hello

numUnique <- hello %>% 
  group_by(Sample, Genus)

library(phyloseq)

# collapse to rank (e.g. "Genus")
ps_rank <- Leray_ps
ps_rank <- tax_glom(Leray_ps, taxrank = "Species")   # replace "Genus" with your rank

# get OTU matrix and ensure taxa are rows
mat <- as(otu_table(ps_rank), "matrix")
if (!taxa_are_rows(ps_rank)) mat <- t(mat)

# count taxa with > 0 reads per sample
counts_per_sample <- colSums(mat > 0)

# put into a data.frame
df_counts <- data.frame(Sample = names(counts_per_sample),
                        n_taxa = as.integer(counts_per_sample),
                        row.names = NULL,
                        stringsAsFactors = FALSE)
df_counts <- df_counts %>% arrange(desc(n_taxa))

length((get_taxa_unique(Leray_ps, taxonomic.rank = "Class")))


tax_df <- as.data.frame(as.matrix(tax_table(Leray_ps)), stringsAsFactors = FALSE)
tax_df$ASV <- taxa_names(Leray_ps) 
tax_df %>%
  filter(!is.na(Species) & Species != "") %>%
  group_by(Species) %>%
  summarise(n_ASVs = n()) %>%
  arrange(desc(n_ASVs))



seqtab <- otu_table(Leray_ps)
table(colSums(seqtab>0))
table(rowSums(seqtab>0))

morph_names <- sample_names(morph_ps)
comp_folmer_ps <- subset_samples(Folmer_ps, Station_ID %in% morph_names)


sample_data(comp_folmer_ps)
comp_folmer_ps <- merge_samples(comp_folmer_ps, "Station_ID")
sample_data(comp_folmer_ps)

sample_names(comp_folmer_ps)
sort(sample_names(morph_ps)) == sort(sample_names(comp_folmer_ps))


#################################################################
# Required packages
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
# optional for a quick pairs plot
if (!requireNamespace("GGally", quietly = TRUE)) install.packages("GGally")
library(GGally)

# 1) helper: count unique taxa per sample at a given rank
count_taxa_per_sample <- function(ps, rank = "Genus", min_abundance = 1,
                                  include_unclassified = FALSE) {
  # ps: phyloseq object
  # rank: string, e.g. "Genus" or "Species"
  # min_abundance: minimum reads to consider taxon present in sample
  # include_unclassified: TRUE treat NA/"" as "Unclassified", FALSE drop them
  df <- psmelt(ps) %>%
    filter(Abundance >= min_abundance) %>%
    mutate(taxon = as.character(.data[[rank]])) %>%
    mutate(taxon = ifelse(is.na(taxon) | taxon == "",
                          ifelse(include_unclassified, "Unclassified", NA_character_),
                          taxon)) %>%
    filter(!is.na(taxon)) %>%
    group_by(Sample) %>%
    summarise(n_taxa = n_distinct(taxon), .groups = "drop")
  
  # ensure every sample in ps is present (0 if none)
  all_samples <- sample_names(ps)
  out <- tibble::tibble(Sample = all_samples) %>%
    left_join(df, by = "Sample") %>%
    mutate(n_taxa = ifelse(is.na(n_taxa), 0L, as.integer(n_taxa)))
  
  out
}

# 2) assemble counts from a named list of phyloseq objects
compute_counts_matrix <- function(ps_list, rank = "Genus", min_abundance = 1,
                                  include_unclassified = FALSE) {
  # ps_list: named list, e.g. list(ps1 = ps1, ps2 = ps2, ps3 = ps3, ps4 = ps4)
  if (is.null(names(ps_list)) || any(names(ps_list) == "")) {
    stop("ps_list must be a named list (names will become column names).")
  }
  
  # take sample set from first object; ensure others contain same sample names
  samples_ref <- sample_names(ps_list[[1]])
  # optional: check that all sample sets match (warning if not)
  for (nm in names(ps_list)[-1]) {
    if (!setequal(samples_ref, sample_names(ps_list[[nm]]))) {
      warning("Sample names differ between first object and ", nm,
              ". Intersection will be used.")
    }
  }
  
  # compute counts per object
  counts_df <- tibble::tibble(Sample = samples_ref)
  for (nm in names(ps_list)) {
    temp <- count_taxa_per_sample(ps_list[[nm]], rank, min_abundance, include_unclassified)
    # align rows by Sample (use intersection)
    temp <- temp %>% filter(Sample %in% samples_ref)
    counts_df <- counts_df %>% left_join(temp, by = "Sample")
    colnames(counts_df)[ncol(counts_df)] <- nm
  }
  # ensure integer columns
  for (nm in names(ps_list)) counts_df[[nm]] <- as.integer(counts_df[[nm]])
  counts_df
}

# 3) Quick pairs plot (GGally)
plot_counts_pairs_ggpairs <- function(counts_df, log_transform = FALSE) {
  # counts_df: output of compute_counts_matrix (Sample column + object columns)
  mat <- counts_df %>% select(-Sample)
  if (log_transform) mat <- log10(mat + 1)
  # GGally will show histograms on diagonal, scatter+smooth and correlations
  p <- GGally::ggpairs(mat,
                       upper = list(continuous = wrap("cor", size = 4)),
                       lower = list(continuous = wrap("smooth", alpha = 0.3, se = FALSE)),
                       diag = list(continuous = "barDiag"))
  p
}

# 4) Custom pairwise ggplot with lm + R^2 annotation (saves list of plots)
pairwise_scatter_plots <- function(counts_df, log_transform = FALSE, prefix = "pairplot") {
  # returns a named list of ggplot objects
  mat <- counts_df %>% select(-Sample)
  if (log_transform) {
    mat <- mat %>% mutate_all(~log10(.x + 1))
  }
  
  obj_names <- colnames(mat)
  combos <- combn(obj_names, 2, simplify = FALSE)
  
  plots <- list()
  for (cmb in combos) {
    a <- cmb[1]; b <- cmb[2]
    dat <- tibble::tibble(x = mat[[a]], y = mat[[b]], Sample = counts_df$Sample)
    fit <- lm(y ~ x, data = dat)
    r2 <- summary(fit)$r.squared
    p <- ggplot(dat, aes(x = x, y = y)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
      labs(x = a, y = b, title = paste(a, "vs", b),
           subtitle = paste0("R² = ", formatC(r2, digits = 3))) +
      theme_minimal()
    name <- paste(a, "_vs_", b, sep = "")
    plots[[name]] <- p
  }
  plots
}

#########################
# Example usage:
# ps_list <- list(obj1 = ps1, obj2 = ps2, obj3 = ps3, obj4 = ps4)
# counts_df <- compute_counts_matrix(ps_list, rank = "Genus", min_abundance = 1, include_unclassified = FALSE)
# View(counts_df)  # Sample | obj1 | obj2 | obj3 | obj4
#
# Quick matrix plot:
# p_pairs <- plot_counts_pairs_ggpairs(counts_df, log_transform = FALSE)
# print(p_pairs)
#
# Or get individual pair plots:
# pair_plots <- pairwise_scatter_plots(counts_df, log_transform = TRUE)
# # show one:
# print(pair_plots[[1]])
# # save all:
# for (nm in names(pair_plots)) ggsave(filename = paste0("plots/", nm, ".png"), plot = pair_plots[[nm]], width = 5, height = 5)
#########################
sample_names(morph_ps)

sample_data(morph_ps)    
sample_data(Leray_ps)      
test <- merge_samples(Leray_ps, group = "Station_ID")
sample_data(test)
