library(dada2); packageVersion("dada2")
library(argparser, quietly=TRUE)

# Usage: dada_tax_assign.R --RDS --refFile --spref --outdir --outprefix

##### SCRIPT MODE #####
p <- arg_parser("DADA2 in R")

# Defining arguments
p = add_argument(p, "--RDS", help="Enter seqtabnochim RDS filepath", type='character')
p = add_argument(p, "--refFile", help="Enter assignTaxonomy ref file", type='character')
p = add_argument(p, "--spref", help="Enter addSpecies ref file", type='character')
p = add_argument(p, "--outdir", help="Enter output directory", type='character')
p = add_argument(p, "--outprefix", help = "Enter output prefix for filenames", type = 'character')

argv <- parse_args(p)

# Extracting argument values (Script mode)
RDS_path = argv$RDS
spref = argv$spref
refDB = argv$refFile
outdir = argv$outdir
outprefix = paste0(outdir, "/", argv$outprefix, "_")

#####

##### INTERACTIVE MODE #####
# RDS_path <-  "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/03_dadaASV/dada_asv_output/mICOintF/mICOintF_seqtab.nochim.Rds"
# spref <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/DBs/final_dada_formatted/id_cluster_noempty_newparse.clustered.fasta"
# refDB <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/DBs/final_dada_formatted/tax_cluster_noempty_newparse.clustered.fasta"
# outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/04_dadaTax/dada_tax_output/mICOintF1"
# outprefix <- "mICOintF"
# outprefix <- paste0(outdir, "/", outprefix, "_")

#####

dir.create(outdir)

# Load in files
print("Loading files")
seqtab.nochim <- readRDS(RDS_path)
otutab <- t(seqtab.nochim)

# Assign taxonomy
print("assignTaxonomy")
taxa <- assignTaxonomy(seqtab.nochim, refDB, multithread = TRUE, tryRC = TRUE)

print("addSpecies")
taxa <- addSpecies(taxa, spref)

# Writing outputs
print("Writing outputs")
write.table(taxa, file = paste0(outprefix, "taxa.tsv"), quote = FALSE)
saveRDS(taxa, paste0(outprefix, "taxa.Rds"))

# Create abundance table with tax assignments and counts - merged by ASVs (rownames)
abundance <- merge(taxa, otutab, by = "row.names", all = TRUE)
write.csv(abundance, file = paste0(outprefix, "abundance.csv"), quote = FALSE)

# The following line replaces any instances of NNNNNNNNNN with nothing
# Use this if paired reads didn't merge in dada2, and you used justConcatenate = TRUE
# abundance$Row.names <- gsub("NNNNNNNNNN", "", abundance$Row.names)
