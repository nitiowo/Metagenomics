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

setwd('/Volumes/Samsung_1TB/Zooplankton/')
outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/alpha_diversity/"
marker <- "18S"
markerdir <- paste0(outdir, marker)
dir.create(markerdir, recursive = TRUE)

####################################################################

# Importing data and generating basic data objects

####################################################################

SampleMetadata <- read.csv('/Volumes/Samsung_1TB/Zooplankton/resources/data/zoop96_metadata.csv', header = TRUE, stringsAsFactors = TRUE)
seqtab.nochim <- readRDS('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_SSUF/seqtab.nochim.Rds')
tax_info <- readRDS('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_SSUF/taxa.Rds')
abundance <- readRDS('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_SSUF/abundance_SSUF.Rds')

# asv_ref <- read.csv('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_LCO/ASV_ref.csv')
# asv_fasta <- read.delim('/Volumes/Samsung_1TB/Zooplankton/zoop_backup/dada_outputs/dada_output_LCO/ASVs.fa')

rownames(SampleMetadata) <- rownames(seqtab.nochim)
SampleMetadata$Mesh <- as.factor(SampleMetadata$Mesh)
SampleMetadata$Year_Sampled <- as.factor(SampleMetadata$Year_Sampled)
SampleMetadata$Tube <- as.factor(SampleMetadata$Tube)

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

# For 18S ----
tax_info1 <- as.matrix(readRDS("/Volumes/Samsung_1TB/Zooplankton/phyloseq/SSUF/wrangling/SSU_animalia_taxFix_taxClean.RDS"))

ps_temp <- prune_taxa(rownames(tax_info1), ps)
ps_temp1 <- phyloseq(otu_table(otu_table(ps_temp)), 
                     tax_table(tax_info1),
                     refseq(refseq(ps_temp)),
                     sample_data(sample_data(ps_temp)))
ps <- ps_temp1

####################################################################

# Plotting rarefaction curve

####################################################################

# We are only plotting the rarefaction curve to visualize ASV read depth
# We are not actually rarefying data

# Run this entire chunk together
png(filename = paste0(markerdir, "/", marker, "_rarecurve.png"), 
    width = 1000, height = 800)
rarecurve(seqtab.nochim, step=100, lwd=2, ylab="ASVs", label=F)
abline(v=(min(rowSums(seqtab.nochim))), col = "gray") # adding a vertical line at the fewest seqs in any sample
title(main = paste0(marker, " rarefaction curve"))
dev.off()

####################################################################

# Creating subset data for comparisons

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

# Collapsing data to Species level. NAs removed by default. Use NArm = FALSE if you want to keep
ps.Species <- tax_glom(ps, taxrank = "Species")



####################################################################

# Setting filename variables for plots

####################################################################

# Use these variables to set the subset being plotted, and change filenames and plot titles accordingly
ps <- ps
level <- "ASV"
tax_subset <- "unfiltered"

# Create subdirectories if they don't exist
outfilepath <- paste0(outdir, marker, "/", level, "/", tax_subset, "/")
dir.create(outfilepath, recursive = TRUE)

# Create output file prefix
outfilename <- paste0(outfilepath, marker, "_", level, "_", tax_subset)

####################################################################

# Plotting alpha diversity boxplots

####################################################################

# Plotting richness and diversity by Lake ----

# Set the comparisons you are interested in here:
a_my_comparisons <- list( c("Erie", "Huron"), 
                          c("Huron", "Michigan"), 
                          c("Ontario", "Superior"), 
                          c("Superior", "Michigan"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                   symbols = c("****", "***", "**", "*", "ns"))


# Actual plot call:
lake_rich.p <- plot_richness(ps, x="Lake", color = "Lake", 
                             measures = c("Observed", "Simpson")) +
  geom_boxplot(alpha=0.6) + theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45, 
                                                         hjust=1, 
                                                         vjust=1, 
                                                         size=12)) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, 
                     label = "p.signif", symnum.args = symnum.args)+
  ggtitle(paste0(marker, " - Alpha diversity measures by Lake - ", 
                 level, " - ", tax_subset)) +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
lake_rich.p

# Save plot
ggsave(filename = paste0(outfilename, "_aDiv_byLake.png"),
       plot = last_plot(), width = 1000, height = 800, device = png, 
       units = "px", scale = 3)


# Plotting richness and diversity by site by Mesh ----
mesh_rich.p <- plot_richness(ps, x="Mesh", measures = c("Observed", "Simpson")) +
  geom_boxplot(alpha=0.6) + theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45, 
                                                         hjust=1, 
                                                         vjust=1, 
                                                         size=12)) +
  ggtitle(paste0(marker, " - Alpha diversity measures by Mesh - ", 
                 level, " - ", tax_subset)) +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
mesh_rich.p

# Save plot
ggsave(filename = paste0(outfilename, "_aDiv_byMesh.png"), 
       plot = last_plot(), width = 1000, height = 800, device = png, 
       units = "px", scale = 3)


# Plotting richness and diversity by site by Benthic ----
benth_rich.p <- plot_richness(ps, x="Morph_avail", measures = c("Observed", "Simpson")) +
  geom_boxplot(alpha=0.6) + theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45, 
                                                         hjust=1, 
                                                         vjust=1, 
                                                         size=12)) +
  ggtitle(paste0(marker, " - Alpha diversity measures by Morophological data availability - ", 
                 level, " - ", tax_subset)) +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
benth_rich.p

# Save plot
ggsave(filename = paste0(outfilename, "_aDiv_byBenth.png"), 
       plot = last_plot(), width = 1000, height = 800, device = png, 
       units = "px", scale = 3)




####################################################################

# Calculating and saving statistical test outputs

####################################################################


# Calculate alpha diversity metrics for each sample
richness <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
write.csv(richness, file = paste0(outfilename, "_aDiv_metrics.csv"))

# Run anova and see if variable (Lake, Mesh, etc) affects diversity metric
anova.sh = aov(richness$Simpson ~ sample_data(ps)$Lake)
summary(anova.sh)

# Based on anova results, we can calculate Tukey SDs (or Kruskal-Wallis test if non-normal data)
TukeyHSD(anova.sh)
kruskal.test(richness$Simpson ~ sample_data(ps)$Lake)

# Get list of p-values from wilcoxon test for each pairwise comparison
pairwise.wilcox.test(richness$Simpson, sample_data(ps)$Lake, p.adj = "bonf")

# Output statistical test summaries to txt file
capture.output(
  summary(anova.sh),
  TukeyHSD(anova.sh),
  kruskal.test(richness$Simpson ~ sample_data(ps)$Lake),
  pairwise.wilcox.test(richness$Simpson, sample_data(ps)$Lake, p.adj = "bonf"), 
  file = paste0(outfilename, "_aDiv_byLake_statTests.txt"))

          