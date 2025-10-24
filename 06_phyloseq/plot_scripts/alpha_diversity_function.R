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
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

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

ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/morph_ps_integer.RDS")

outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/alpha_diversity_function_test2/"
marker <- "Morphological" # Only works with no slash at the end
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

# Alpha diversity boxplot function ----

####################################################################

plot_adivBox <- function(ps, group_variables, comparisons, measures,
                     marker_outdir, marker, level, tax_subset){
  
  for(group_var in group_variables){
    
    if(group_var == "Morph_avail" && grepl("Benthic", tax_subset)){
      next
    } else if (marker == "Morphological" && group_var != "Lake"){
      next
    }
    
    print(paste0("Plotting alpha diversity ", marker, " at ", level, " level by ", group_var,".png"))
    
    outfilepath <- paste0(marker_outdir, "/", level, "/", tax_subset, "/")
    outfilepath1 <- paste0(outfilepath,"/adiv_box/")
    dir.create(outfilepath1, recursive = TRUE)
    outfilename <- paste0(outfilepath1, marker, "_", level, "_", tax_subset)
    
    # Significance thresholds
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                       symbols = c("****", "***", "**", "*", "ns"))
    
    # Actual plot call:
    rich.p <- plot_richness(ps, x=group_var, color = group_var, 
                                 measures = measures) +
      geom_boxplot(alpha=1, outlier.shape = NA) + theme_bw() + 
      geom_jitter(width = 0.2, height = 0.2, alpha = 0.4) +
      theme(legend.position="none", 
            axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12)) +
      stat_compare_means(method = "wilcox.test", comparisons = comparisons, 
                         label = "p.signif", symnum.args = symnum.args) +
      ggtitle(paste0(marker, " - Alpha diversity measures by ", group_var," - ", 
                     level, " - ", tax_subset)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Show plot
    rich.p$layers <- rich.p$layers[-1]
    rich.p
    
    # Save plot
    ggsave(filename = paste0(outfilename, "_aDiv_by", group_var, ".png"),
           plot = rich.p, width = 1200, height = 800, device = png, 
           units = "px", scale = 3)
  }
}

# Arguments for interactive function run

# a_my_comparisons <- list( c("Erie", "Huron"), 
#                           c("Huron", "Michigan"), 
#                           c("Ontario", "Superior"), 
#                           c("Superior", "Michigan"))
# adiv_measures <- c("Observed", "Shannon", "InvSimpson")
# group_vars <- c("Lake", "Mesh", "Morph_avail")
# level <- "ASV"
# tax_subset <- "unfiltered"


plot_adivBox(ps, group_vars, a_my_comparisons, adiv_measures,
             markerdir, marker, level, tax_subset)


####################################################################

# Global Indicators ---- 

####################################################################

# Create subdirectories if they don't exist

global_adiv <- function(ps,
                        marker_outdir, marker, level, tax_subset){
  
  outfilepath <- paste0(marker_outdir, "/", level, "/", tax_subset, "/")
  outfilepath1 <- paste0(outfilepath, "global_aDiv/")
  dir.create(outfilepath1, recursive = TRUE)
  
  # Create output file prefix
  outfilename <- paste0(outfilepath1, marker, "_", level, "_", tax_subset)
  
  print(paste0("Writing ", basename(outfilename), "_aDiv_metrics.csv"))
  
  tab <- microbiome::alpha(ps, index = "all")
  
  write.csv(tab, file = paste0(outfilename, "_global-aDiv-indicators.csv"), 
            col.names = T, row.names = T)
  
  return(tab)
}
  
global_adiv(ps, markerdir, marker, level, tax_subset)


####################################################################

# Alpha diversity statistical tests ---- 

####################################################################

adiv_stats <- function(ps, group_variables, measures,
                       marker_outdir, marker, level, tax_subset){
  
  outfilepath <- paste0(marker_outdir, "/", level, "/", tax_subset, "/")
  outfilepath1 <- paste0(outfilepath, "aDiv_stats/")
  dir.create(outfilepath1, recursive = TRUE)
  
  for(group_var in group_variables){
    
    if(group_var == "Morph_avail" && grepl("Benthic", tax_subset)){
      next
    } else if (marker == "Morphological" && group_var != "Lake"){
      next
    }
    
    # Create output file prefix
    outfilename <- paste0(outfilepath1, marker, "_", level, "_", tax_subset)
    
    print(paste0("Writing ", basename(outfilename), "_aDiv_stats_by", group_var, ".txt"))
    
    # Calculate alpha diversity metrics for each sample
    richness <- estimate_richness(ps, measures = measures)
    
    i <- 0
    for(measure in measures){
      
      # Run anova and see if variable (group_var, Mesh, etc) affects diversity metric
      anova.sh = aov(richness[[measure]] ~ sample_data(ps)[[group_var]])
      summary(anova.sh)
      
      # Based on anova results, we can calculate Tukey SDs (or Kruskal-Wallis test if non-normal data)
      TukeyHSD(anova.sh)
      kruskal.test(richness[[measure]] ~ sample_data(ps)[[group_var]])
      
      # Get list of p-values from wilcoxon test for each pairwise comparison
      pairwise.wilcox.test(richness[[measure]], sample_data(ps)[[group_var]], p.adj = "bonf")
      
      # Overwrite file on first loop, then append thereafter
      if(i == 0){
        boo <- FALSE
      } else{
        boo <- TRUE
      }
      
      # Output statistical test summaries to txt file
      capture.output(
        cat(paste0(">ANOVA summary for ", measure, " by ", group_var, ":")),
        summary(anova.sh),
        
        cat(paste0("\n>Tukey HSD based on ANOVA for ", measure, " by ", group_var, ":")),
        TukeyHSD(anova.sh),
        
        cat(paste0("\n>Kruskal test for ", measure, " by ", group_var, ":")),
        kruskal.test(richness[[measure]] ~ sample_data(ps)[[group_var]]),
        
        cat(paste0("\n>Pairwise Wilcoxon test for ", measure, " by ", group_var, ":")),
        pairwise.wilcox.test(richness[[measure]], sample_data(ps)[[group_var]], p.adj = "bonf"), 
        
        file = paste0(outfilename, "_aDiv_stats_by", group_var, ".txt"),
        append = boo)
      
      i <- i + 1
    }
  }
}

adiv_stats(ps, group_vars, measures,
           markerdir, marker, level, tax_subset)

####################################################################

# Creating master plot function ----

####################################################################

# Select this whole code chunk before running it
plot_adiv <- function(ps_object, marker, tax_levels, subsets, marker_outdir,
                      box_array, stats_array){
  # Removing problem samples
  if (marker == "Folmer"){
    ps1 <- subset_samples(ps_object, !Tube %in% c("68","9"))
  } else if (marker == "18S"){
    ps1 <- subset_samples(ps_object, Tube != "57")
  } else if (marker == "Leray" || marker == "Morphological"){
    ps1 <- ps_object
  }
  
  # Iterating through taxa levels ----
  ####################################################################
  for(level in tax_levels){
    
    print(paste0("Processing at ", level, " level"))
    
    # Collapsing data to specified tax_level. NAs removed by default. Use NArm = FALSE if you want to keep
    if(level != "ASV"){
      ps2 <- tax_glom(ps1, taxrank = level)
    } else if(level == "ASV"){
      ps2 <- ps1
    }
    
    # Iterating through subsets  ----
    ####################################################################
    for(tax_subset in subsets){
      
      print(paste0("Processing at ", level, " level for ", tax_subset, " subset"))
      
      if(marker == "Morphological" && grepl("Benthic", tax_subset)){
        next
      }
      
      ####################################################################
      
      # Creating subset data for comparisons ----
      
      ####################################################################
      
      if(tax_subset == "noRotifer"){
        # Subsetting out Rotifers
        ps <- subset_taxa(ps2, Phylum != "Rotifera")
      } else if(tax_subset == "onlyRotifer"){
        # Subsetting only Rotifers
        ps <- subset_taxa(ps2, Phylum == "Rotifera")
      } else if(tax_subset == "noBenthic"){
        # Subsetting samples with Morphological data
        ps <- subset_samples(ps2, Morph_avail == "YES")
      } else if(tax_subset == "onlyBenthic"){
        # Subsetting samples without Morphological data
        ps <- subset_samples(ps2, Morph_avail == "NO")
      } else if(tax_subset == "noRotifer_noBenthic"){
        # Subsetting out samples with no Morphological data and Rotifers
        ps <- ps2 %>% 
          subset_taxa(Phylum != "Rotifera") %>% 
          subset_samples(Morph_avail == "YES")
      } else if(tax_subset == "unfiltered"){
        ps <- ps2
      }
      
      # Create subdirectories if they don't exist
      outfilepath <- paste0(marker_outdir, "/", level, "/", tax_subset, "/")
      dir.create(outfilepath, recursive = TRUE)
      
      # Create output file prefix
      outfilename <- paste0(outfilepath, marker, "_", level, "_", tax_subset)
      
      # Sanity check
      print(paste("ntaxa = ", ntaxa(ps), "nsamples = ", nsamples(ps)))
      
      ##########################################
      # Alpha diversity  boxplots ----
      
      groupsy <- box_array$Group_vars
      comps <- box_array$Comparisons
      measures <- box_array$Measures
      
      plot_adivBox(ps, groupsy, comps, measures, 
               marker_outdir, marker, level, tax_subset)
      
      ##########################################
      # Global metrics ----
      
      global_adiv(ps, 
                  marker_outdir, marker, level, tax_subset)
      
      ##########################################
      # Alpha Diversity statistical tests ----
      
      groupsy <- stats_array$Group_vars
      measures <- stats_array$Measures
      
      adiv_stats(ps, groupsy, measures,
                 markerdir, marker, level, tax_subset)
    }
  }
}


####################################################################

# Calling function (creating plots) ----

####################################################################
tax_levels <- c("ASV", "Species", "Genus")
subsets <- c("unfiltered", "noRotifer", "noBenthic", "noRotifer_noBenthic")

# Available subsets: unfiltered, noRotifer, noBenthic, noRotifer_noBenthic, onlyRotifer, onlyBenthic
# Available group variables: Lake, Mesh, Morph_avail
# Available measures: ("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")


my_comparisons <- list( c("Erie", "Huron"),
                          c("Huron", "Michigan"),
                          c("Ontario", "Superior"),
                          c("Superior", "Michigan"))

box_array <- list(Comparisons = my_comparisons, 
                       Group_vars = c("Lake", "Mesh", "Morph_avail"), 
                       Measures = c("Observed", "Shannon", "Simpson", "Chao1"))

stats_array <- list(Group_vars = c("Lake", "Mesh", "Morph_avail"), 
                  Measures = c("Observed", "Shannon", "Simpson", "Chao1"))


# Function call
plot_adiv(ps, marker, tax_levels, subsets, markerdir,
          box_array, stats_array)

