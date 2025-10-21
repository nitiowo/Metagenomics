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
outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/beta_diversity/"


####################################################################

# Importing data and generating basic data objects 

####################################################################

SampleMetadata <- read.csv('/Volumes/Samsung_1TB/Zooplankton/resources/data/zoop96_metadata.csv', header = TRUE, stringsAsFactors = TRUE)
seqtab.nochim <- readRDS("/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_LCO/seqtab.nochim.Rds")
tax_info <- as.matrix(readRDS("/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_LCO/taxa.Rds"))

# abundance <- readRDS('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_SSUF/abundance_SSUF.Rds')
# asv_ref <- read.csv('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_LCO/ASV_ref.csv')
# asv_fasta <- read.delim('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_LCO/ASVs.fa')

rownames(SampleMetadata) <- rownames(seqtab.nochim)
SampleMetadata$Mesh <- as.factor(SampleMetadata$Mesh)
SampleMetadata$Year_Sampled <- as.factor(SampleMetadata$Year_Sampled)
SampleMetadata$Tube <- as.factor(SampleMetadata$Tube)

# Creating phyloseq object
SAM <- sample_data(SampleMetadata)
OTU <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
TAX <- tax_table(tax_info)

ps <- phyloseq(SAM, OTU, TAX)

sample_names(ps) <- paste0("S-", sample_names(ps))

# Rename ASV sequences to headers
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
ps

# Or import ps object directly ----
ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/morph_ps_integer.RDS")

marker <- "Morphological"
markerdir <- paste0(outdir, marker)
dir.create(markerdir, recursive = TRUE)

####################################################################

# Creating subset data for comparisons ----

####################################################################

# Making subsets
if (marker == "Folmer"){
  ps <- subset_samples(ps, !Tube %in% c("68","9"))
} else if (marker == "18S"){
  ps <- subset_samples(ps, Tube != "57")
}

# Subsetting Rotifers
ps.noRotifer <- subset_taxa(ps, Phylum != "Rotifera")
ps.Rotifer <- subset_taxa(ps, Phylum == "Rotifera")

# Subsetting samples with Morphological data
ps.noBenthic <- subset_samples(ps, Morph_avail == "YES")
ps.Benthic <- subset_samples(ps, Morph_avail == "NO")

# Subsetting out samples with no Morphological data and Rotifers
ps.noRotifer.noBenthic <- ps %>% 
  subset_taxa(Phylum != "Rotifera") %>% 
  subset_samples(Morph_avail == "YES")

# Collapsing data to Species level ----
# NAs removed by default. Use NArm = FALSE if you want to keep
ps.Species <- tax_glom(ps, taxrank = "Species")
ps_temp <- ps
ps <- ps.Species

# Subsetting Rotifers
ps.noRotifer.Species <- subset_taxa(ps, Phylum != "Rotifera")
ps.Rotifer.Species <- subset_taxa(ps, Phylum == "Rotifera")

# Subsetting samples with Morphological data
ps.noBenthic.Species <- subset_samples(ps, Morph_avail == "YES")
ps.Benthic.Species <- subset_samples(ps, Morph_avail == "NO")

# Subsetting out samples with no Morphological data and Rotifers
ps.noRotifer.noBenthic.Species <- ps %>% 
  subset_taxa(Phylum != "Rotifera") %>% 
  subset_samples(Morph_avail == "YES")

ps <- ps_temp

# ps_list <- list(c(ps,
#                   ps.noBenthic, ps.noRotifer, ps.noRotifer.noBenthic,
#                   ps.Species, ps.noRotifer.Species, ps.noBenthic.Species, ps.noRotifer.noBenthic.Species))
# 
# ps_hash <- hash()
# ps_hash[["ps"]] <- list(level = "ASV", subset = "unfiltered")
# ps_hash[["ps"]] <- list(level = "ASV", subset = "unfiltered")
# ps_hash[["ps"]] <- list(level = "ASV", subset = "unfiltered")
# ps_hash[["ps"]] <- list(level = "ASV", subset = "unfiltered")
# ps_hash[["ps"]] <- list(level = "ASV", subset = "unfiltered")
# ps_hash[["ps"]] <- list(level = "ASV", subset = "unfiltered")
# ps_hash[["ps"]] <- list(level = "ASV", subset = "unfiltered")
# ps_hash[["ps"]] <- list(level = "ASV", subset = "unfiltered")

####################################################################

# Setting filename variables for plots ----

####################################################################

# Use these variables to set the subset being plotted, and change filenames and plot titles accordingly
ps <- ps.noRotifer
level <- "Species"
tax_subset <- "noRotifer"

# Create subdirectories if they don't exist
outfilepath <- paste0(outdir, marker, "/", level, "/", tax_subset, "/")
dir.create(outfilepath, recursive = TRUE)

# Create output file prefix
outfilename <- paste0(outfilepath, marker, "_", level, "_", tax_subset)

####################################################################

# NMDS plot 

####################################################################

# Transform data to proportions as appropriate for Bray-Curtis distances (NMDS)
# See for info: https://jkzorz.github.io/2019/06/06/NMDS.html
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))  # Transforms sample counts to relative abundance. Not sure why we do this
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")  # NMDS needs an ordination object to plot. Default distance measure is bray-curtis

# Plot NMDS
lake_nmds.p <- plot_ordination(ps.prop, ord.nmds.bray, color="Lake") + 
  stat_ellipse(aes(group=Lake)) +
  geom_text(aes(label=sample_names(ps), hjust=0.3, vjust=-0.4), size = 3) +
  ggtitle(paste0(marker, " - NMDS ordination by Lake - ", 
                 level, " - ", tax_subset))

# Show plot
lake_nmds.p

# Save plot
ggsave(filename = paste0(outfilename, "_NMDS_byLake.png"),
       plot = last_plot(), width = 1000, height = 800, device = png, 
       units = "px", scale = 3)

####################################################################

# PCoA plot 

####################################################################

# Ordinate and plot using PCoA
ord = ordinate(ps.prop, method="PCoA", distance = "bray")

# Plot PCoA
lake_pcoa.p <- plot_ordination(ps, ord, color = "Lake") + 
  geom_point(size=2) + 
  stat_ellipse(aes(group=Lake)) +
  geom_text(aes(label=sample_names(ps), hjust=0.3, vjust=-0.4)) +
  ggtitle(paste0(marker, " - PCoA ordination by Lake - ", 
                 level, " - ", tax_subset)) +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
lake_pcoa.p

# Save plot
ggsave(filename = paste0(outfilename, "_PCoA_byLake.png"),
       plot = last_plot(), width = 1000, height = 800, device = png, 
       units = "px", scale = 3)

####################################################################

# Plot phylogenetic tree

####################################################################

# Phylogenetic tree
random_tree <- rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

dev.off()
png(filename = paste0(outdir, "/randtree.png"), width = 1000, height = 800)
plot(random_tree, color = "Lake")
dev.off()

ps1 <- merge_phyloseq(ps, random_tree) # Add tree to ps object

####################################################################

# Heatmap and Dendrogram

####################################################################

# Heatmap
plot_heatmap(ps1)

# Save plot
ggsave(filename = paste0(outfilename, "_heatmap.png"),
       plot = last_plot(), width = 1000, height = 800, device = png, 
       units = "px", scale = 3)

# Dendrogram
asv_tab <- as.data.frame(t(otu_table(ps)))

# For morphological: asv_tab <- as.data.frame(otu_table(ps))

samp <- as.data.frame(sample_data(ps))

deseq_counts <- DESeqDataSetFromMatrix(asv_tab, colData = samp, design = ~Lake) 

deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
transOTU <- otu_table(vst_trans_count_tab, taxa_are_rows = TRUE)

sample_names(transOTU) <- sample_names(ps)
ps.DEtrans <- merge_phyloseq(transOTU, sample_data(ps), tax_table(ps))

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
#plot(euc_clust) 00B0F6

#show_col(hue_pal()(5))
colorCodes <- c(Erie="#F8766D", Superior="#E76BF3", Huron="#A3A500", Ontario="#00B0F6", Michigan = "#00BF7D")

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
labels_colors(euc_dend) <- colorCodes[samp$Lake][order.dendrogram(euc_dend)]

# Run this whole chunk together
png(filename = paste0(outfilename, "_dendrogram.png"), width = 800, height = 800)
plot(euc_dend)
title(ylab="Euclidean distance", 
      main = paste0(marker, " - Dendrogram by Lake - ", level, " - ", tax_subset), 
      cex.main = 1.9, cex.axis = 2)
legend("topright", legend=c("Erie", "Superior", "Huron", "Ontario", "Michigan"),
      col=colorCodes, pch=19, cex = 1.5, pt.cex = 1, y.intersp = 1, x.intersp = 1)
dev.off()



