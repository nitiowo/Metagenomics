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
library(RColorBrewer)

setwd('/Volumes/Samsung_1TB/Zooplankton/')
outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/beta_diversity_function_test/"

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

orig_ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/morph_ps_integer.RDS")
ps <- orig_ps

outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/beta_diversity_function_test/"

marker <- "Morphological" # Only works with no slash at the end
markerdir <- paste0(outdir, marker)
dir.create(markerdir, recursive = TRUE)


random_tree <- rtree(ntaxa(orig_ps), rooted=TRUE, 
                     tip.label=taxa_names(orig_ps))
ps <- merge_phyloseq(orig_ps, random_tree) # Add tree to ps object


####################################################################

# Ordination plots 

####################################################################

ord_plot <- function(ps, ord_methods, distance_metrics, group_variables,
                     marker_outdir, marker, level, tax_subset){
  
  for(meth in ord_methods){
    for(dist in distance_metrics){
      for(group_var in group_variables){
        
        if(group_var == "Morph_avail" && grepl("Benthic", tax_subset)){
          next
        } else if (marker == "Morphological" && group_var != "Lake"){
          next
        }
        
        print(paste("Plotting", meth, dist, group_var, sep = " "))
        
        outfilepath <- paste0(marker_outdir, "/", level, "/", tax_subset, "/")
        
        outfilepath1 <- paste0(outfilepath, meth, "/", dist, "/")
        dir.create(outfilepath1, recursive = TRUE)
        outfilename <- paste0(outfilepath1, marker, "_", level, "_", tax_subset)
        
        if(dist == "jaccard"){
          # Jaccard uses only presence/absence data (not abundances). 
          # Transformation to binary not strictly necessary (ordinate function natively uses only binary data for Jaccard)
          ps <- transform_sample_counts(ps, function(x) ifelse(x > 0, 1L, 0L))
          ord <- ordinate(ps, method=meth, distance="jaccard")
        } else if(dist == "bray"){
          # Transform data to relative abundance as appropriate for Bray-Curtis distances
          # See for info: https://jkzorz.github.io/2019/06/06/NMDS.html
          ps <- transform_sample_counts(ps, function(otu) otu/sum(otu))
          ord <- ordinate(ps, method=meth, distance="bray")  
        }
        
        # Plot Ordination
        ord.p <- plot_ordination(ps, ord, color=group_var) + 
          stat_ellipse(aes(group=.data[[group_var]])) +
          geom_text(aes(label=sample_names(ps), hjust=0.3, vjust=-0.4), size = 3) +
          ggtitle(paste0(marker, " - ", dist, " ", meth, " ordination by ", 
                         group_var," - ", level, " - ", tax_subset))
        
        # Show plot
        ord.p
        
        # Save plot
        ggsave(filename = paste0(outfilename, "_", dist,"_", meth, "_by", group_var,".png"),
               plot = ord.p, width = 1000, height = 800, device = png, 
               units = "px", scale = 3)
        
      }
    }
  }
}


dists <- c("bray", "jaccard")
groupsy <- c("Lake", "Mesh")
ords <- c("PCoA", "NMDS")
level <- "Species"
tax_subset <- "unfiltered"

ord_plot(ps, ords, dists, groupsy, 
         markerdir, marker, level, tax_subset)


####################################################################

# Heatmap

####################################################################

heat_plot <- function(ps, distance_metrics, group_variables,
                      marker_outdir, marker, level, tax_subset){
  
    for(dist in distance_metrics){
      for(group_var in group_variables){
        
        if(group_var == "Morph_avail" && grepl("Benthic", tax_subset)){
          next
        } else if (marker == "Morphological" && group_var != "Lake"){
          next
        }
        
        print(paste("Plotting heatmap", dist, group_var, sep = " "))
        
        outfilepath <- paste0(marker_outdir, "/", level, "/", tax_subset, "/")
        
        outfilepath1 <- paste0(outfilepath,"/heatmap/", dist, "/")
        dir.create(outfilepath1, recursive = TRUE)
        outfilename <- paste0(outfilepath1, marker, "_", level, "_", tax_subset)
        
        if(level == "ASV"){
          heatlevel <- NULL 
        } else {
            heatlevel <- level
          }
        
        heat.p <- plot_heatmap(ps, taxa.label = heatlevel, 
                               sample.order = group_var, 
                          distance = dist) +
          ggtitle(paste0(marker, " - ", dist, " heatmap by ", 
                         group_var," - ", level, " - ", tax_subset))
        
        ggsave(filename = paste0(outfilename, "_", dist, "_by", group_var,"_heatmap.png"),
               plot = heat.p, width = 1000, height = 800, device = png, 
               units = "px", scale = 3)
      }
    }
}


dists <- c("jaccard", "bray")
groupsy <- c("Lake")
ords <- c("PCoA", "NMDS")
level <- "ASV"
tax_subset <- "noRotifer"

heat_plot(ps, dists, groupsy, 
         markerdir, marker, level, tax_subset)


####################################################################

# Dendrogram

####################################################################

dend_plot <- function(ps, group_variables,
                      marker_outdir, marker, level, tax_subset){
  
  asv_tab <- as.data.frame(t(otu_table(ps)))
  
  #For morphological: 
  if(marker == "Morphological"){
    asv_tab <- as.data.frame(otu_table(ps))
  }
  
  samp <- as.data.frame(sample_data(ps))
  for(group_var in group_variables){
    
    if(group_var == "Morph_avail" && grepl("Benthic", tax_subset)){
      next
    } else if (marker == "Morphological" && group_var != "Lake"){
      next
    }
    
    print(paste("Plotting dendrogram", group_var, sep = " "))
    
    outfilepath <- paste0(marker_outdir, "/", level, "/", tax_subset, "/")
    
    outfilepath1 <- paste0(outfilepath,"/dendrogram/")
    dir.create(outfilepath1, recursive = TRUE)
    outfilename <- paste0(outfilepath1, marker, "_", level, "_", tax_subset)
    
    deseq_counts <- DESeqDataSetFromMatrix(asv_tab, colData = samp, design = reformulate(group_var)) 
    
    deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
    deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
    
    # and here is pulling out our transformed table
    vst_trans_count_tab <- assay(deseq_counts_vst)
    
    # and calculating our Euclidean distance matrix
    euc_dist <- dist(t(vst_trans_count_tab))
    
    euc_clust <- hclust(euc_dist, method="ward.D2")
    
    # hclust objects like this can be plotted with the generic plot() function
    #plot(euc_clust) 00B0F6
    
    #show_col(hue_pal()(5))
    
    make_colors <- function(groupvar, pal = "Set2"){
      lv <- unique(as.character(groupvar))
      lv[is.na(lv)] <- "NA"                 # treat NA as its own level (optional)
      n <- length(lv)
      maxc <- brewer.pal.info[pal, "maxcolors"]
      cols <- if(n <= maxc) brewer.pal(n, pal) else colorRampPalette(brewer.pal(maxc, pal))(n)
      setNames(cols, lv)
    }
    
    colorCodes <- make_colors(samp[[group_var]], pal="Set2")
    legendnames <- unique(as.character(samp[[group_var]]))
    
    euc_dend <- as.dendrogram(euc_clust, hang=0.1)
    labels_colors(euc_dend) <- colorCodes[samp[[group_var]]][order.dendrogram(euc_dend)]
    
    # Run this whole chunk together
    png(filename = paste0(outfilename, "_by", group_var,"_dendrogram.png"), width = 800, height = 800)
    plot(euc_dend)
    title(ylab="Euclidean distance", 
          main = paste0(marker, " - Dendrogram by ", group_var, " - ", level, " - ", tax_subset), 
          cex.main = 1.9, cex.axis = 2)
    legend("topright", legend=legendnames, 
           col=colorCodes, pch=19, cex = 1.5, pt.cex = 1, y.intersp = 1, x.intersp = 1)
    dev.off()
  }
}


groupsy <- c("Lake", "Mesh")
level <- "ASV"
tax_subset <- "unfiltered"

dend_plot(ps, groupsy, 
          marker_outdir, marker, level, tax_subset)


####################################################################

# Creating plot function ----

####################################################################

# Select this whole code chunk before running it
plot_bdiv <- function(ps_object, marker, tax_levels, subsets, marker_outdir,
                      ord_plot_array, heat_plot_array, dend_plot_array){
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
      # Ordination Plots ----
      
      dists <- ord_plot_array$Distances
      groupsy <- ord_plot_array$Group_vars
      ords <- ord_plot_array$Ord_Methods
      
      ord_plot(ps, ords, dists, groupsy, 
               marker_outdir, marker, level, tax_subset)
      
      ##########################################
      # Heatmap Plots ----
      
      dists <- heat_plot_array$Distances
      groupsy <- heat_plot_array$Group_vars
      
      heat_plot(ps, dists, groupsy, 
                marker_outdir, marker, level, tax_subset)
      
      ##########################################
      # Dendrogram Plots ----
      
      groupsy <- dend_plot_array$Group_vars
      
      dend_plot(ps, groupsy, 
                marker_outdir, marker, level, tax_subset)
    }
  }
}


####################################################################

# Calling function (creating plots) ----

####################################################################
tax_levels <- c("ASV", "Species", "Genus")
subsets <- c("unfiltered", "noRotifer", "noBenthic", "noRotifer_noBenthic")

# Available subsets: unfiltered, noRotifer, noBenthic, noRotifer_noBenthic, onlyRotifer, onlyBenthic
# Available distances: jaccard, bray
# Available group variables: Lake, Mesh, Morph_avail
# Available ordination methods: PCoA, NMDS

ord_plot_array <- list(Distances = c("jaccard", "bray"), 
                       Group_vars = c("Lake", "Mesh", "Morph_avail"), 
                       Ord_Methods = c("PCoA", "NMDS"))

heat_plot_array <- list(Distances = c("jaccard", "bray"), 
                       Group_vars = c("Lake", "Mesh", "Morph_avail"))

dend_plot_array <- list(Group_vars = c("Lake", "Mesh", "Morph_avail"))

# Calling plot function
plot_bdiv(ps, marker, tax_levels, subsets, markerdir,
          ord_plot_array, heat_plot_array, dend_plot_array)
