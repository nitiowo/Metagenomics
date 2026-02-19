library(dada2); packageVersion("dada2")
library(argparser, quietly=TRUE)

# Usage: dada_tax_assign.R --RDS --refFile --outdir --outprefix --runspecies [--spref]

##### SCRIPT MODE #####
p <- arg_parser("DADA2 in R")

# Defining arguments
p = add_argument(p, "--RDS", help="Enter seqtabnochim RDS filepath", type='character')
p = add_argument(p, "--refFile", help="Enter assignTaxonomy ref file", type='character')
p = add_argument(p, "--spref", help="Enter addSpecies ref file (required if --runspecies true)", type='character', default=NA)
p = add_argument(p, "--outdir", help="Enter output directory", type='character')
p = add_argument(p, "--outprefix", help = "Enter output prefix for filenames", type = 'character')
p = add_argument(p, "--runspecies", help = "Run addSpecies? (true/false)", type = 'character', default = "true")

argv <- parse_args(p)

# Extracting argument values (Script mode)
RDS_path = argv$RDS
refDB = argv$refFile
outdir = argv$outdir
outprefix = paste0(outdir, "/", argv$outprefix, "_")
run_species = tolower(argv$runspecies) == "true"
spref = argv$spref

#####

##### INTERACTIVE MODE #####
# RDS_path <-  "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/03_dadaASV/dada_asv_output/mICOintF/mICOintF_seqtab.nochim.Rds"
# spref <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/DBs/final_dada_formatted/id_cluster_noempty_newparse.clustered.fasta"
# refDB <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/DBs/final_dada_formatted/tax_cluster_noempty_newparse.clustered.fasta"
# outdir <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/04_dadaTax/dada_tax_output/mICOintF1"
# outprefix <- "mICOintF"
# outprefix <- paste0(outdir, "/", outprefix, "_")
# run_species <- TRUE
#####

dir.create(outdir)

# Load in files
print("Loading files")
seqtab.nochim <- readRDS(RDS_path)
otutab <- t(seqtab.nochim)

# Assign taxonomy
print("assignTaxonomy")
taxa <- assignTaxonomy(seqtab.nochim, refDB, multithread = TRUE, tryRC = TRUE)

# Optionally run addSpecies if --runspecies is true
if (run_species) {
  if (is.na(spref)) stop("--spref is required when --runspecies is true")
  print("addSpecies")
  taxa <- addSpecies(taxa, spref)
} else {
  print("Skipping addSpecies (runspecies = false)")
}

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
