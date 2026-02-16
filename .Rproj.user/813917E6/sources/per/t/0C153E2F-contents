library(readr)
library(dplyr)
library(stringr)
library(utils)
library(ggplot2)
library(forcats)
library(taxize)

setwd('/Volumes/Samsung_1TB/Zooplankton/')

total_bold <- read.csv("./Metagenomics/DBs/BOLD_pull_direct/COIonly_taxa_table.tsv", fill = TRUE, header = F, sep = "\t") 
outdir <- "./Metagenomics/DBs/DB_wrangling/"
treb <- read.csv("./Metagenomics/DBs/DB_wrangling/Trebitz_list_2025.txt", fill = TRUE, header = TRUE, sep = "\t") #Trebitz taxa list

treb.unique.Class <- unique(treb$Class)
treb.unique.Order <- unique(treb$Order)
treb.unique.Species <- unique(treb$Species)

# Reduce to unique rows, add count column
result <- total_bold %>%
  select(-1) %>%                     # drop ID column (first column)
  count(across(everything(), name = "count") %>%
  arrange(desc(count))

# Change col names
cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species", "Subspecies", "count")
colnames(result) <- cols
cols1 <- c("ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species", "Subspecies")
colnames(total_bold) <- cols1

# Grab only Trebitz classes from the BOLD database
relevant_classes <- unique(c("Copepoda", "Thecostraca", "Malacostraca", "Tantulocarida", 
                      "Maxillopoda", "Monogononta", "Ichthyostraca", treb.unique.Class))

only_relevant_from_BOLD <- result %>% 
  filter(Class %in% treb.unique.Class |
           Class %in% relevant_classes)

# Grab sequence IDs from original BOLD dataset using taxa in only_relevant dataset
withID_only_relevant <- total_bold %>% 
  filter(Class %in% treb.unique.Class |
           Class %in% relevant_classes)

# Grab rows with weird species
bad_rows <- withID_only_relevant %>%
  filter(grepl("[^A-Za-z ]", Species))

# Grab weird rows with full species info
# unicode-aware pattern: letters (any language), apostrophe and hyphen allowed
pat <- "^[\\p{L}'-]+\\s+[\\p{L}'-]+$"

ID_relevant_subspeciesSaved <- withID_only_relevant %>%
  # ensure character types so if_else doesn't choke on factors
  mutate(
    Species = as.character(Species),
    Subspecies = if ("Subspecies" %in% names(.) ) as.character(Subspecies) else NA_character_
  ) %>%
  mutate(
    trimmed = str_trim(Species),
    token_count = str_count(trimmed, "\\S+"),
    first_two = str_extract(trimmed, "^\\S+\\s+\\S+"),
    # everything after the first two tokens (trim and convert empty -> NA)
    extra_info = str_trim(str_replace(trimmed, "^\\S+\\s+\\S+", "")),
    extra_info = na_if(extra_info, ""),
    
    # flag rows for which we want to perform the replacement (>=3 tokens and first_two are letters)
    flag = token_count >= 3 & str_detect(first_two, regex(pat, ignore_case = TRUE)),
    
    # perform replacements only for flagged rows
    Subspecies = if_else(flag, extra_info, Subspecies),
    Species    = if_else(flag, first_two, Species)
  ) %>%
  # remove helper columns
  select(-trimmed, -token_count, -first_two, -flag, -extra_info, -Subfamily)

write.csv(ID_relevant_subspeciesSaved, paste0(outdir, "relevant_taxa_IDs.csv"))

##### Random sampling of outgroups ----
# Agglomerate to specific level, remove relevant classes
Phylum_glom_count <- total_bold %>% 
  filter(!Class %in% relevant_classes) %>% 
  group_by(Phylum) %>% 
  summarise(Count = n()) %>% 
  arrange(desc(Count))

# Random sampling for top 10 classes
classes <- head(Class_glom_count$Class, n = 10)

set.seed(42)
sampled_IDs <- total_bold %>% 
  filter(Class %in% classes) %>% 
  group_by(Class) %>% 
  slice_sample(n = 100) %>% 
  ungroup()

write.csv(sampled, paste0(outdir, "outgroups_top10Classes.csv"))


# Species level IDs
species_level_IDs <- ID_relevant_subspeciesSaved %>% 
  filter(Species != "None" & !is.na(Species)) %>% 
  select(ID, Species)

species_sampled <- sampled %>% 
  filter(Species != "None" & !is.na(Species)) %>% 
  select(ID, Species)

species_withOutgroups <- rbind(species_level_IDs, species_sampled)
write.csv(species_withOutgroups, paste0(outdir, "speciesIDs_withOutgroups.csv"))




non_unique <- (species_withOutgroups[duplicated(species_withOutgroups),])
count(non_unique)

non_unique1 <- species_level_IDs[duplicated(species_level_IDs),]
non_unique2 <- species_sampled[duplicated(species_sampled),]

non_unique3 <- withID_only_relevant[duplicated(withID_only_relevant), ]







