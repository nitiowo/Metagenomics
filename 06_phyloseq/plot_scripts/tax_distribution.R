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

orig_ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/Leray_dadaOriginal.RDS")
ps <- orig_ps

outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/tax_distribution/"

marker <- "Leray" # Only works with no slash at the end
markerdir <- paste0(outdir, marker)
dir.create(markerdir, recursive = TRUE)


random_tree <- rtree(ntaxa(orig_ps), rooted=TRUE, 
                     tip.label=taxa_names(orig_ps))
ps <- merge_phyloseq(orig_ps, random_tree) # Add tree to ps object



####################################################################

# Taxonomic distribution

####################################################################

tax_plot <- function(ps, top_x, group_variables,levels,
                     marker_outdir, marker, tax_subset){
  for(tax_subset in tax_subsets){
    for(top in top_x){
      for(level in levels){
        for(group_var in group_variables){
          
          if(group_var == "Morph_avail" && grepl("Benthic", tax_subset)){
            next
          } else if (marker == "Morphological" && group_var != "Lake"){
            next
          }
          
          topchar <- as.character(top)
          
          print(paste("Plotting top", top, level, group_var, sep = " "))
          
          outfilepath <- paste0(marker_outdir, "/", tax_subset, "/")
          
          outfilepath1 <- paste0(outfilepath, level, "/top", topchar, "/")
          dir.create(outfilepath1, recursive = TRUE)
          
          outfilebase <- paste0(outfilepath1, marker, "_", level, "_", tax_subset)
          outfilename <- paste0(outfilebase, "_", level, "_top", top, "_by", group_var,".png")
          
          
          title <- paste0(marker, " Taxonomic distribution of top ", top, ' ', 
                          level, " by ", group_var, " - ", tax_subset)
          # Plot bar chart of top 20 abundant ASVs
          topx <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
          ps.transform <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
          ps.topx <- prune_taxa(topx, ps.transform)
          
          topx_glom <- tax_glom(ps.topx, taxrank = level)
          
          bar.p <- plot_bar(topx_glom, fill=level, x="Station_ID") + 
            #x = "Station_ID") +
            # theme_bw() +
            {if(group_var != "all")facet_wrap(vars(.data[[group_var]]), scales="free_x")} +
          #facet_wrap(~Lake) + 
          ggtitle(title)  +
            theme(plot.title = element_text(size = 25, hjust = 0.5), 
                  axis.text.x = element_text(angle = 90)) +
            ylab("Relative abundance") 
          # geom_bar(stat="identity", position = "stack")
          
          # scale_fill_discrete(breaks = c("Huron", "Michigan", "Superior", "Erie", "Ontario")) # change order of legend
          show(bar.p)
          return(bar.p)
          
          # plot_bar(ps.transform, fill="Order") +
          #   facet_wrap(~Lake, scales = "free_x") +
          #   ylab("Relative abundance") + 
          #   ggtitle("18S - Overall biodiversity of Great Lakes at Order level from all sequences")
          
          # ggsave using last_plot() to save whatever versions you want. Make sure to change names
          # ggsave(filename = paste0(outdir, "/morph_station_top20_wt.png"), plot = last_plot(), 
          #        width = 1000, height = 800, device = png, units = "px", scale = 3)
          
        }
      }
    }
  }
}

# top_xs <- c(10, 20, 40, "all")
# groups <- c("Lake", "Mesh", "Station_ID")
# ranks <- c("Species", "Genus", "Family", "Order", "Class")
# tax_subsets <- c("unfiltered", "noRotifer", "noBenthic", "noRotifer_noBenthic")

top_xs <- c(10)
groups <- c("Mesh")
ranks <- c("Family")
tax_subsets <- c("unfiltered")


tax_plot(ps, top_x = top_xs, group_variables = groups, levels = ranks,
        marker_outdir = markerdir, marker = marker, tax_subset = tax_subsets)
  

  
  
  
  
paste0("haha", top)
  
tax_subset <- "unfiltered"
group_var <- "Lake"
level <- "Species"
top <- 20
topchar <- as.character(top)
  
outfilepath <- paste0(markerdir, "/", "unfiltered", "/")

outfilepath1 <- paste0(outfilepath, "Species", "/top", topchar, "/")
dir.create(outfilepath1, recursive = TRUE)

outfilebase <- paste0(outfilepath1, marker, "_", level, "_", tax_subset)
outfilename <- paste0(outfilebase, "_", level, "_top", top, "_by", group_var,".png")
  
  
  
  
  
  
  
  
  
  


