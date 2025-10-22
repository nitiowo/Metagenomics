library(dada2)
library(stringr)
library(Biostrings)
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd('/Volumes/Samsung_1TB/Zooplankton/Metagenomics/DBs/')

# read fasta
seqs <- readDNAStringSet("bold_clustered.assigntaxonomy.fasta")

# Identify headers without "_"
has_underscore <- grepl("_", names(seqs))
no_underscore  <- !has_underscore

# Subset sequences
seqs_with_underscore <- seqs[has_underscore]
seqs_without_underscore <- seqs[no_underscore]

# Subset sequences
seqs_filtered <- seqs[keep]
length(seqs)           # original
length(seqs_without_underscore)  # filtered
length(seqs_with_underscore) #filtered out

# extract headers
headers <- names(seqs)

# Split headers into taxonomic ranks
tax_list <- str_split(headers, ";")

# # Determine max rank depth
# max_depth <- max(sapply(tax_list, length))
# max_depth <- max_depth-1
# 
# # Pad shorter entries with NA so all have same length
# tax_matrix <- t(sapply(tax_list, function(x) {
#   length(x) <- max_depth-1
#   x
# }))

tax_matrix <- do.call(rbind, tax_list)

# Optional: convert to data frame
tax_df <- as.data.frame(tax_matrix, stringsAsFactors = FALSE)
# colnames(tax_df) <- paste0("Rank", 1:7)
tax_df <- tax_df[1:(length(tax_df)-1)]

# Name columns
colnames(tax_df) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Convert to matrix
#tax_matrix <- as.matrix(tax_matrix)

# # Use Species column (or other rank)
# species_names <- tax_matrix[, "Species"]
# 
# # Keep unique species only
# unique_species <- unique(species_names)
# unique_rows <- unique(tax_matrix)

# # Count how many sequences per taxa
# species_counts <- table(tax_matrix)
# count_df <- data.frame(species_counts)

tax_df <- tax_df %>%
  group_by(across(everything())) %>%   # group by all 7 columns
  summarise(count = n(), .groups="drop") %>%
  arrange(desc(count))

# Make a simple OTU table for phyloseq
# OTU <- otu_table(matrix(species_counts, ncol=1, dimnames=list(unique_species, "Reference")), taxa_are_rows=TRUE)
# TAX <- tax_table(tax_matrix[match(unique_species, tax_matrix[, "Species"]), , drop=FALSE])
# phy_ref <- phyloseq(OTU, TAX)
# 
# 
# # Create an OTU table with counts = 1 for each sequence (just for phyloseq structure)
# otu <- matrix(1, nrow=length(seqs), ncol=1)
# rownames(otu) <- names(seqs)
# colnames(otu) <- "Reference"
# 
# # Create phyloseq objects
# OTU <- otu_table(otu, taxa_are_rows=TRUE)
# TAX <- tax_table(tax_matrix)
# 
# phy_ref <- phyloseq(OTU, TAX)

# Function to plot taxonomic distribution
plot_tax_distribution <- function(tax_df, rank = "Rank6", top_n = 20, rotate_x = 45, title = NULL) {
  # Check rank exists
  if (!rank %in% colnames(tax_df)) {
    stop(paste("Rank column", rank, "not found in tax_df"))
  }
  if (!"count" %in% colnames(tax_df)) {
    stop("tax_df must have a 'count' column")
  }
  
  # Summarize counts at the chosen rank
  tax_summary <- tax_df %>%
    group_by(.data[[rank]]) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    arrange(desc(total_count)) %>%
    slice_max(total_count, n = top_n)  # keep top N
  
  # Default title if none provided
  if (is.null(title)) {
    title <- paste("Top", top_n, "taxa at", rank)
  }
  
  # Create bar plot
  p <- ggplot(tax_summary, aes(x = reorder(.data[[rank]], -total_count), 
                               y = total_count, fill = .data[[rank]])) +
    geom_col() +
    labs(x = rank, y = "Number of reference sequences", title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = rotate_x, hjust = 1),
          legend.position = "none")
  
  return(p)
}



# tax_df = your data frame with Rank1-Rank7 + count
plot_tax_distribution(tax_df, rank = "Species", top_n = 50)  # Phylum
plot_tax_distribution(tax_df_with, rank = "Species", top_n = 50)  # Genus
plot_tax_distribution(tax_df_without, rank = "Class", top_n = 50)  # Species

nrow(tax_df_with[tax_df_with$Phylum=="Rotifera",])

noinsect_without <- tax_df_without[tax_df_without$Class!="Insecta",]
plot_tax_distribution(noinsect_without, rank = "Class", top_n = 50)  # Species

uniq_spec <- unique(tax_df$Species)

#############
# TAXIZE #
#############


result.long <- tax_df$Species %>%
  gnr_resolve(data_source_ids = c(1,3,8,9,11,12), 
              with_canonical_ranks=T)

d1 <- split(unique_genus_species, ceiling(seq_along(unique_genus_species)/1000))

result.long <- vector(mode = "list", length = length(d1))
for (i in 1:length(d1)){
  result.long[[i]] <- d1[[i]] %>%
    gnr_resolve(data_source_ids = c(1,3,8,9,11,12), 
                with_canonical_ranks=T)
  print(i)
}

result.long <- data.table::rbindlist(result.long)

tst <- as.data.frame(d1[[1]])

result.short <- result.long %>%
  select(submitted_name, matched_name2, score)%>%
  distinct()
result.short


write.table(result.short,
            "./taxize_outputs/result.short.BOLD.txt", 
            sep="\t", row.names = F, quote = F)

# Go through the above text file, add 3 columns: implement, alternative, dupl. 
# Default: all implement = TRUE, all alternative = empty, all dupl = FALSE
# Save in new file with .COMMENTS.

corr.df.BOLD <- read.table("./taxize_outputs/result.short.BOLD.COMMENTS.txt", 
                           sep="\t", header=T, stringsAsFactors = F)
corr.df.BOLD %<>% filter(!dupl ==T)


# Turning species names to lowercase to fix case change issues.
# Creating duplicate species name column too keep originals

corr.df.BOLD$dupname <- corr.df.BOLD$submitted_name
comb_mICOint$dupSpec <- comb_mICOint$Species

comb_mICOint <- comb_mICOint %>% mutate(Species = tolower(Species))
corr.df.BOLD <- corr.df.BOLD %>% mutate(matched_name2 = tolower(matched_name2))


# Join new names to original abundance table
mICOint_join <- comb_mICOint %>% 
  # join operation:
  left_join(corr.df.mICOint, by=c("Species" ="submitted_name")) %>%
  mutate(new.taxon = ifelse(implement == T, matched_name2,
                            ifelse(implement == F & is.na(alternative)==T, taxa,
                                   ifelse(implement == F & is.na(alternative)==F, alternative))))

write.table(taxa.df.2,
            "tax.fixed.mICOint.txt", 
            sep="\t", row.names = F, quote = F)


# Reading fixed abundance tables back in for use in analysis
COI_join <- read.table("tax.fixed.COI.txt", header = TRUE, sep = "\t")

# Subsetting only taxa names
COI.tax <- pull(COI_join, new.taxon)
BOLD.tax <- pull(corr.df.BOLD, matched_name2)


saveRDS(BOLD.tax, file = './taxize_outputs/BOLD.tax.Rds')










