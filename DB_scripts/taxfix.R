library(taxize)
library(readr)
library(dplyr)
library(utils)
library(ggplot2)
library(forcats)

setwd('/Volumes/Samsung_1TB/Zooplankton/')

total_bold <- read.csv("./Metagenomics/DBs/BOLD_pull_direct/taxa_table.tsv", fill = TRUE, header = F, sep = "\t") 
outdir <- "./Metagenomics/DBs/DB_wrangling/"
treb <- read.csv("./Metagenomics/DBs/DB_wrangling/Trebitz_list_2025.txt", fill = TRUE, header = TRUE, sep = "\t") #Trebitz taxa list

##### Generating agglomerated taxa tables ----
treb.unique.Class <- unique(treb$Class)
treb.unique.Order <- unique(treb$Order)
treb.unique.Species <- unique(treb$Species)

result <- total_bold %>%
  select(-1) %>%                     # drop ID column (first column)
  count(across(everything()), name = "count") %>%
  arrange(desc(count))

cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species", "Subspecies", "count")
colnames(result) <- cols

only_relevant_from_BOLD <- result %>% 
  filter(Class %in% treb.unique.Class |
           Class %in% c("Copepoda", "Thecostraca", "Malacostraca", "Tantulocarida", "Maxillopoda", "Monogononta", "Ichthyostraca"))


classy_fixed_onlyrelevant <- data.frame(matrix(ncol = 7, nrow = 0))
ranks <- c("Kingdom", "Phylum", "Class", "Order", 
           "Family", "Genus", "Species")
colnames(classy_fixed_onlyrelevant) <- ranks
uniq_genus <- unique(only_relevant_from_BOLD$Genus)  
# ----

##### Function definitions ----
collapse_to_rank_unique <- function(df, rank,
                                    ranks_order = c("Kingdom","Phylum","Class","Order","Family",
                                                    "Subfamily", "Genus","Species","Subspecies")) {
  # find which of the canonical ranks are present in df (preserve canonical order)
  present <- ranks_order[tolower(ranks_order) %in% tolower(names(df))]
  # make sure requested rank is present
  if(!tolower(rank) %in% tolower(present)) stop("Requested rank not found in dataframe")
  # keep columns up to and including the requested rank
  keep <- present[1:which(tolower(present) == tolower(rank))]
  # map to original column names (preserve exact case from df)
  keep <- names(df)[match(tolower(keep), tolower(names(df)))]
  df %>% select(all_of(keep)) %>% distinct()
}

collapsed <- collapse_to_rank_unique(only_relevant_from_BOLD, "Genus")

relevant_genus_glom <- collapsed %>% 
  select(-Subfamily)
  
fix_class <- function(tax_info.df, start_rank){
  
  # Making empty dataframe to store results
  new_tax <- as.data.frame(matrix(NA, nrow = nrow(tax_info.df), ncol = ncol(tax_info.df)))
  row.names(new_tax) <- c(rownames(tax_info.df))
  colnames(new_tax) <- colnames(tax_info.df)
  
  # Create empty vector to store highest assigned value
  new_tax.list <- vector(mode = "list", length = nrow(tax_info.df))
  
  # Initialize number of ncoll number to iterate through (8 because 8 levels of taxonomy)
  badcalls <- matrix(nrow = nrow(tax_info.df)*8, ncol = 3)
  colnames(badcalls) <- c("asv_numbers", "taxa", "Fixed?")
  
  # Iterate through each row (each taxon)
  for(i in 1:nrow(tax_info.df)){ 
    curr_asv <- rownames(tax_info.df)[i]
    ncoll <- as.numeric(which(colnames(tax_info.df) == start_rank))
    
    # Iterate through each column (each rank)
    for(j in 1:ncoll){
      tax <- tax_info.df[i,ncoll]  # Start in start_rank column, save taxonomic name to variable
      curr_rank <- colnames(tax_info.df[ncoll])
      
      if(is.na(tax) || grepl("unculture|environment|None", tax)){  # if tax is NA or contains "uncultured" or "environmental", go back one column and check Genus, etc...
        ncoll <- ncoll - 1  # If it's NA, change column to 6 and check Genus, etc...
        print(paste0(curr_asv," Going back because ", curr_rank, " NA"))
        next
        
        # Attempting classification
      } else {
        print(paste0(curr_asv, " Attempting classification for ", curr_rank, " ", tax))
        
        db_info <- (get_gbifid_(tax))[[1]] # Change get_x_ function to preferred db.
        
        if(nrow(db_info) == 0){
          ncoll <- ncoll - 1
          print(paste0("No classificaiton results for ", tax))
          badcalls[i+(j-1), 1] <- curr_asv  # Saving the ASV headers and tax names. 
          badcalls[i+(j-1), 2] <- tax # If there is 
          next
        }
        db_rank <- db_info %>% 
          subset(rank == tolower(curr_rank) &
                   status == "ACCEPTED" &
                   (matchtype == "EXACT" | matchtype == "HIGHERRANK"))
        db_name <- db_rank$canonicalname
        
        tryCatch({
          classy <- classification(db_name, db = "worms")
          print(classy)# When you find one that's not NA, attempt classification
          classy <- as.data.frame(matrix(unlist(classy), ncol = 3))

          if(is.na(classy$V1[1])){
            print("No WORMS result. Trying ITIS")
            classy <- classification(db_name, db = "itis")   # When you find one that's not NA, attempt classification
            classy <- as.data.frame(matrix(unlist(classy), ncol = 3))  # Unlist results to access variables
          }
          
          if(is.na(classy[tolower(classy$V2) == "kingdom", 2]) || # If it's NA at Kingdom level
             (nrow(classy[tolower(classy$V2) == "kingdom",]) == 0)){  # OR there is no Kingdom row, the taxonomic info was not good enough, so we go back one column
            print(paste0(curr_asv, " Kingdom-level failure - saving to badcalls. Retrying with higher rank"))
            badcalls[i+(j-1), 1] <- curr_asv  # Saving the ASV headers and tax names. 
            badcalls[i+(j-1), 2] <- tax # If there is more than one bad call per row, it saves the second bad call in the next line
            ncoll <- ncoll - 1
            
          } else {   # If taxonomic info is not NA, and it classifies properly, we save the results in new_tax at the levels we want.
            print("Saving values")
            new_tax[i,7] <- ifelse(nrow(classy[tolower(classy$V2) == "species",]) == 0, NA, (classy[tolower(classy$V2) == "species",1]))
            new_tax[i,6] <- ifelse(nrow(classy[tolower(classy$V2) =="genus",]) == 0, NA, (classy[tolower(classy$V2) =="genus",1]))
            new_tax[i,5] <- ifelse(nrow(classy[tolower(classy$V2) =="family",]) == 0, NA, (classy[tolower(classy$V2) =="family",1]))
            new_tax[i,4] <- ifelse(nrow(classy[tolower(classy$V2) =="order",]) == 0, NA, (classy[tolower(classy$V2) =="order",1]))
            new_tax[i,3] <- ifelse(nrow(classy[tolower(classy$V2) =="class",]) == 0, NA, (classy[tolower(classy$V2) =="class",1]))
            new_tax[i,2] <- ifelse(nrow(classy[tolower(classy$V2) =="phylum",]) == 0, NA, (classy[tolower(classy$V2) =="phylum",1]))
            new_tax[i,1] <- ifelse(nrow(classy[tolower(classy$V2) =="kingdom",]) == 0, NA, (classy[tolower(classy$V2) =="kingdom",1]))
            ncoll <- 7
            break
          }
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
  badcalls <- as.data.frame(badcalls)
  badcalls_clean <- badcalls[rowSums(is.na(badcalls)) != ncol(badcalls), ]
  return(list(new_tax = new_tax, badcalls = badcalls_clean))
}
# ----

subset <- relevant_genus_glom[18:30,]

obj <- fix_class(relevant_genus_glom, "Genus")
newtax <- obj$new_tax
badcalls <- obj$badcalls

write.csv(newtax, paste0(outdir, "newtax.csv"))
write.csv(badcalls, paste0(outdir, "badcalls.csv"))


classification("Limnadopsis", db = "worms")


# classy <- classification("ouwghwohig", db = "worms") 
# str(classy)
# print(classy)# When you find one that's not NA, attempt classification
# classy <- as.data.frame(matrix(unlist(classy), ncol = 3)) 
# 
# classy
# 
# is.na(classy$V1[1])


