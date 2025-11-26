library(readr)
library(dplyr)
library(utils)
library(ggplot2)
library(forcats)
library(taxize)

setwd('/Volumes/Samsung_1TB/Zooplankton/')

total_bold <- read.csv("./Metagenomics/DBs/BOLD_pull_direct/taxa_table.tsv", fill = TRUE, header = F, sep = "\t") 
outdir <- "./Metagenomics/DBs/DB_wrangling//"
treb <- read.csv("./Metagenomics/DBs/DB_wrangling/Trebitz_list_2025.txt", fill = TRUE, header = TRUE, sep = "\t") #Trebitz taxa list

treb.unique.Class <- unique(treb$Class)
treb.unique.Order <- unique(treb$Order)
treb.unique.Species <- unique(treb$Species)


# df <- read_tsv("taxa_table.tsv", col_types = cols(.default = "c"))  # all cols as character

result <- total_bold %>%
  select(-1) %>%                     # drop ID column (first column)
  count(across(everything()), name = "count") %>%
  arrange(desc(count))

cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species", "Subspecies", "count")
colnames(result) <- cols
# 
# result %>%
#   group_by(Phylum) %>%
#   summarise(total = sum(Count, na.rm = TRUE)) %>%
#   arrange(desc(total)) %>%
#   mutate(Phylum = forcats::fct_reorder(Phylum, total)) %>%
#   ggplot(aes(x = Phylum, y = total)) +
#   geom_col(fill = "steelblue") +
#   coord_flip() +
#   labs(x = "Phylum", y = "Number of sequences",
#        title = "Distribution of sequences by Phylum") +
#   theme_minimal(base_size = 14)



# Read the collapsed table (expects a column named 'count')
#taxa <- read_tsv("taxa_table_collapsed.tsv", col_types = cols(.default = "c")) %>%
 # mutate(count = as.numeric(count))

# Function: plot subrank (e.g., Class) within rank (e.g., Phylum)
plot_subrank_within_rank <- function(df,
                                     rank = "Phylum",
                                     subrank = "Class",
                                     top_n_rank = 10,    # keep top X phyla
                                     top_n_sub = 6,      # keep top Y classes per phylum
                                     other_label = "Other",
                                     title = NULL,
                                     save_plot = NULL) {
  # defensive: replace NA/empty with label
  df <- df %>%
    mutate(across(all_of(c(rank, subrank)), ~ ifelse(is.na(.) | . == "" , "Unclassified", .)))
  
  # 1) choose top N groups at the 'rank' level by summed counts
  top_ranks <- df %>%
    group_by(across(all_of(rank))) %>%
    summarise(rank_total = sum(count, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(rank_total)) %>%
    slice_head(n = top_n_rank) %>%
    pull(!!sym(rank))
  
  df_top <- df %>% filter(.data[[rank]] %in% top_ranks)
  
  # 2) for each chosen rank, keep top_n_sub subrank values, collapse others to 'Other'
  df_collapsed <- df_top %>%
    group_by(across(all_of(c(rank, subrank)))) %>%
    summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
    group_by(across(all_of(rank))) %>%
    mutate(rank_sub_rank = row_number(desc(total)),
           keep = rank_sub_rank <= top_n_sub,
           !!subrank := ifelse(keep, .data[[subrank]], other_label)) %>%
    ungroup() %>%
    group_by(across(all_of(c(rank, subrank)))) %>%
    summarise(total = sum(total), .groups = "drop")
  
  # 3) order ranks by total for plotting
  rank_order <- df_collapsed %>%
    group_by(across(all_of(rank))) %>%
    summarise(rank_total = sum(total), .groups = "drop") %>%
    arrange(desc(rank_total)) %>%
    pull(!!sym(rank))
  
  df_collapsed <- df_collapsed %>%
    mutate(
      !!rank := factor(.data[[rank]], levels = rank_order),
      !!subrank := fct_infreq(.data[[subrank]])  # subrank factor by overall freq (keeps legend tidy)
    )
  
  # 4) plot
  plt <- ggplot(df_collapsed, aes(x = !!sym(rank), y = total, fill = !!sym(subrank))) +
    geom_col() +
    coord_flip() +
    labs(x = rank, y = "Number of sequences",
         title = ifelse(is.null(title),
                        paste0("Distribution of ", subrank, " within top ", top_n_rank, " ", rank, "s"),
                        title),
         fill = subrank) +
    theme_minimal(base_size = 13)
  
  print(plt)
  
  if (!is.null(save_plot)) {
    ggsave(save_plot, plt, width = 10, height = 6, dpi = 300)
    message("Saved plot to: ", save_plot)
  }
  
  invisible(df_collapsed)  # return the plotting table if user wants to inspect
}


# Collapse identical taxonomy, add count collumn
result <- only_relevant_from_BOLD %>%
  select(-1) %>%                     # drop ID column (first column)
  count(across(everything()), name = "count") %>%
  arrange(desc(count))

# Plot top 12 phyla, keeping top 8 classes per phylum
plot_subrank_within_rank(result,
                         rank = "Class",
                         subrank = "Genus",
                         top_n_rank = 10,
                         top_n_sub = 15,
                         save_plot = paste0(outdir, "phylum_class_stacked.png"))


cols1 <- c("ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species", "Subspecies")
colnames(total_bold) <- cols1

# Subset bold data based on trebitz data
setdiff(treb.unique.Class, unique(total_bold$Class))

# freyi_orders <- total_bold %>% 
#   subset(Species == "Bosmina freyi") %>% 
#   select(Order)
# 
# treb_bdelloids <- treb %>% 
#   subset(Order == "Bdelloidea") %>% 
#   pull(Species)
# 
# setdiff(treb_bdelloids, unique(total_bold$Species))
# 
# bdell_orders <- total_bold %>% 
#   subset(Species == "Dissotrocha aculeata") %>% 
#   pull(Order)


only_relevant_from_BOLD <- total_bold %>% 
  filter(Class %in% treb.unique.Class |
         Class %in% c("Copepoda", "Thecostraca", "Malacostraca", "Tantulocarida", "Maxillopoda", "Monogononta"))

write_csv(only_relevant_from_BOLD, 
          file = paste0(outdir, "BOLD_TrebitzClassSubset_taxList.csv"), col_names = TRUE)

non_unique_elements <- only_relevant_from_BOLD$ID[duplicated(only_relevant_from_BOLD$ID) | duplicated(only_relevant_from_BOLD$ID, fromLast = TRUE)]


length(unique(total_bold$Order))
uniq.treb.species <- (unique(treb$Species))
unique_bold_treb_species <- (unique(only_relevant_from_BOLD$Species))

noBold_treb <- treb$Species[!grepl("v4", treb$BOLD.Entry)] 
Bold_treb <- treb$Species[grepl("v4", treb$BOLD.Entry)] 

missing <- setdiff(uniq.treb.species, unique_bold_treb_species)
missing # What is in the trebitz list, but not the bold subset - missing stuff

no.bold.data <- Reduce(intersect, list(missing, noBold_treb))
no.bold.data     # What's in the list of missing stuff, but also has no bold entry

# Check what's happening with the ones that say they don't have bold entry:
really.notIn.bold <- setdiff(no.bold.data, unique(total_bold$Species))
really.notIn.bold  # Say they don't have BOLD entry, don't appear in master BOLD database.
# These check out - no entry in master BOLD database, or the BOLD website

shouldHave.bold.data <- Reduce(intersect, list(missing, Bold_treb))
shouldHave.bold.data   # Appear in Trebitz, say they have BOLD data, but do not appear in subset.

# Checking if any of the above appear int the master BOLD dataset
fake.notIn.bold <- Reduce(intersect, list(shouldHave.bold.data, unique(total_bold$Species)))
fake.notIn.bold # These appear to be missing from subset and say they have BOLD data according to Trebitz
# But these appear in the master BOLD dataset - so something is stopping them from making it through the filter
# Usually a mismatch in higher level tax assignment (Maxillopoda vs Ichthyostraca)
# Solved by adding relevant classes to the filtering step : Class %in% c("Copepoda", "Thecostraca", "Malacostraca", "Tantulocarida", "Maxillopoda", "Monogononta"))
# We can standardize taxonomy later


# Everything else from shouldHave.bold.data:
fakerm.shouldHave.bold <- setdiff(shouldHave.bold.data, fake.notIn.bold)
fakerm.shouldHave.bold # These should have BOLD data but something weird is happening - mispellings, bad formatting, etc
# If any of these end up having barcodes, we have them in the database.
# Only reason to come back to these is to disambiguate to report how many barcodes were missing
# End goal for this set is to sort it into either 1) Has barcode, and disambiguate name or 2) Doesn't have barcode.

# Sanity check
length(shouldHave.bold.data) == length(fakerm.shouldHave.bold) + length(fake.notIn.bold)









diff2 <- setdiff(diff, uniq.treb.species)
diff2

# result1 <- only_relevant_from_BOLD %>%
#   select(-1) %>%                     # drop ID column (first column)
#   count(across(everything()), name = "count") %>%
#   arrange(desc(count))



"Heteropsyllus nunni" %in% unique_bold_treb_species
"Brachionus bidentata" %in% total_bold$Species


# total_bold$Class[unique(total_bold$Genus)]

noBold_treb <- treb$Species[!grepl("v4", treb$BOLD.Entry)] 

Bold_treb <- treb$Species[grepl("v4", treb$BOLD.Entry)] 

unique(treb$BOLD.Entry[!grepl("v4", treb$BOLD.Entry)])


"Argulus americanus" %in% treb$Species
unique(treb$Class[treb$Species == "Cupelopagis vorax"]) 
#%in% total_bold$Order


"Dicranophorus caudatus" %in% total_bold$Species
unique(total_bold$Class[total_bold$Species == "Cupelopagis vorax"]) 




treb.unique.Class

treb[treb$Class, ]

gna_verifier("Moinidae")

foo <- taxize::get_gbifid("Moinidae")


fool <- classification("Monogononta", db = "gbif")
fool



