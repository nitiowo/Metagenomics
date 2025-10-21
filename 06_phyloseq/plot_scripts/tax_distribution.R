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

# Taxonomic distribution

####################################################################

# Plot bar chart of top 20 abundant ASVs
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.transform <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.transform)

# Plot 
topgen <- subset_taxa(ps.transform, Genus %in% c( "Keratella", 
                                                  "Leptodiaptomus", 
                                                  "Pontellidae_X", 
                                                  "Daphnia", 
                                                  "Synchaeta", 
                                                  "Kellicottia"))

plot_bar(ps.top20, fill="Family", x = "Station_ID") +
  theme_bw() +
  # facet_wrap(~Lake, scales="free_x") 
  #facet_wrap(~Lake) + 
  ggtitle("Morphological")  +
  theme(plot.title = element_text(size = 25, hjust = 0.5), 
        axis.text.x = element_text(angle = 90)) +
  ylab("Relative abundance") +
  geom_bar(stat="identity", position = "stack")

# scale_fill_discrete(breaks = c("Huron", "Michigan", "Superior", "Erie", "Ontario")) # change order of legend


plot_bar(ps.transform, fill="Order") +
  facet_wrap(~Lake, scales = "free_x") +
  ylab("Relative abundance") + 
  ggtitle("18S - Overall biodiversity of Great Lakes at Order level from all sequences")

# ggsave using last_plot() to save whatever versions you want. Make sure to change names
ggsave(filename = paste0(outdir, "/morph_station_top20_wt.png"), plot = last_plot(), 
       width = 1000, height = 800, device = png, units = "px", scale = 3)

