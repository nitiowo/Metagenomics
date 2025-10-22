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

ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/Folmer_dadaOriginal.RDS")

outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/alpha_diversity_function_test/"
marker <- "Folmer" # Only works with no slash at the end
markerdir <- paste0(outdir, marker)
dir.create(markerdir, recursive = TRUE)

####################################################################

# Plotting rarefaction curve ----

####################################################################

# We are only plotting the rarefaction curve to visualize ASV read depth
# We are not actually rarefying data

# Need raw unfiltered seqtab.nochim (asv table) for this. Doing it on subsetted samples is useless.

# Run this entire chunk together
png(filename = paste0(markerdir, "/", marker, "_rarecurve.png"), 
    width = 1000, height = 800)
rarecurve(seqtab.nochim, step=100, lwd=2, ylab="ASVs", label=F)
abline(v=(min(rowSums(seqtab.nochim))), col = "gray") # adding a vertical line at the fewest seqs in any sample
title(main = paste0(marker, " rarefaction curve"))
dev.off()


####################################################################

# Creating plot function ----

####################################################################

plot_adiv <- function(ps_object, marker, tax_levels, subsets, marker_outdir){
  # Removing problem samples
  if (marker == "Folmer"){
    ps <- subset_samples(ps_object, !Tube %in% c("68","9"))
  } else if (marker == "18S"){
    ps <- subset_samples(ps_object, Tube != "57")
  }
  
  # Iterating through taxa levels ----
  ####################################################################
  for(level in tax_levels){
    
    print(paste0("Processing at ", level, " level"))
    
    # Iterating through subsets  ----
    ####################################################################
    for(tax_subset in subsets){
      
      # Collapsing data to specified tax_level. NAs removed by default. Use NArm = FALSE if you want to keep
      if(level != "ASV"){
        ps <- tax_glom(ps_object, taxrank = level)
      } else if(level == "ASV"){
        ps <- ps_object
      }
      
      print(paste0("Processing at ", level, " level for ", tax_subset, " subset"))
      
      ####################################################################
      
      # Creating subset data for comparisons ----
      
      ####################################################################
      
      if(tax_subset == "noRotifer"){
        # Subsetting out Rotifers
        ps <- subset_taxa(ps, Phylum != "Rotifera")
      } else if(tax_subset == "onlyRotifer"){
        # Subsetting only Rotifers
        ps <- subset_taxa(ps, Phylum == "Rotifera")
      } else if(tax_subset == "noBenthic"){
        # Subsetting samples with Morphological data
        ps <- subset_samples(ps, Morph_avail == "YES")
      } else if(tax_subset == "onlyBenthic"){
        # Subsetting samples without Morphological data
        ps <- subset_samples(ps, Morph_avail == "NO")
      } else if(tax_subset == "noRotifer_noBenthic"){
        # Subsetting out samples with no Morphological data and Rotifers
        ps <- ps %>% 
          subset_taxa(Phylum != "Rotifera") %>% 
          subset_samples(Morph_avail == "YES")
      }
      
      # Create subdirectories if they don't exist
      outfilepath <- paste0(marker_outdir, "/", level, "/", tax_subset, "/")
      dir.create(outfilepath, recursive = TRUE)
      
      # Create output file prefix
      outfilename <- paste0(outfilepath, marker, "_", level, "_", tax_subset)
      
      ####################################################################

      # Plotting alpha diversity boxplots ----
      
      ####################################################################
      
      # Plotting richness and diversity by Lake ----
      
      print(paste0("Plotting ", basename(outfilename), "_aDiv_byLake.png"))
      
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
        geom_boxplot(alpha=1, outlier.shape = NA) + theme_bw() + 
        geom_jitter(width = 0.2, height = 0.2, alpha = 0.4) +
        theme(legend.position="none", 
              axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12)) +
        stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, 
                           label = "p.signif", symnum.args = symnum.args) +
        ggtitle(paste0(marker, " - Alpha diversity measures by Lake - ", 
                       level, " - ", tax_subset)) +
        theme(plot.title = element_text(hjust = 0.5))
      
      # Show plot
      lake_rich.p$layers <- lake_rich.p$layers[-1]
      lake_rich.p
      
      # Save plot
      ggsave(filename = paste0(outfilename, "_aDiv_byLake.png"),
             plot = last_plot(), width = 1000, height = 800, device = png, 
             units = "px", scale = 3)
  
      
      # Plotting richness and diversity by site by Mesh ----
      
      if(marker != "Morphological"){
        
        print(paste0("Plotting ", basename(outfilename), "_aDiv_byMesh.png"))
        
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
      }
      
      
      # Plotting richness and diversity by site by Benthic ----
      if(marker != "Morphological" &&
         tax_subset != "noBenthic" && 
         tax_subset != "noRotifer_noBenthic" &&
         tax_subset != "onlyBenthic"){
        
        print(paste0("Plotting ", basename(outfilename), "_aDiv_byBenth.png"))
        
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
      }
      
      
      ####################################################################
      
      # Calculating and saving statistical test outputs ---- 
      
      ####################################################################
      
      print(paste0("Writing ", basename(outfilename), "_aDiv_metrics.csv"))
      
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
      
      print(paste0("Writing ", basename(outfilename), "_aDiv_byLake_statTests.txt"))
      
      # Output statistical test summaries to txt file
      capture.output(
        summary(anova.sh),
        TukeyHSD(anova.sh),
        kruskal.test(richness$Simpson ~ sample_data(ps)$Lake),
        pairwise.wilcox.test(richness$Simpson, sample_data(ps)$Lake, p.adj = "bonf"), 
        file = paste0(outfilename, "_aDiv_byLake_statTests.txt"))
    }
  }
}

####################################################################

# Calling function (creating plots) ----

####################################################################
tax_levels <- c("Species", "Genus")
subsets <- c("noRotifer", "unfiltered", "onlyRotifer")

plot_adiv(ps, marker, tax_levels, subsets, markerdir)

