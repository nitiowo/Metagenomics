library(dada2); packageVersion("dada2")
library(argparser, quietly=TRUE)

p <- arg_parser("DADA2 in R")

p = add_argument(p, "--RDS", help="Enter seqtabnochim RDS filepath", type='character')
p = add_argument(p, "--refFile", help="Enter assignTaxonomy ref file", type='character')
p = add_argument(p, "--spref", help="Enter addSpecies ref file", type='character')
p = add_argument(p, "--outdir", help="Enter output directory", type='character')
p = add_argument(p, "--outprefix", help = "Enter output prefix for filenames", type = 'character')


argv <- parse_args(p)

RDS_path = argv$RDS
spref = argv$spref
refDB = argv$refFile
outdir = argv$outdir
outprefix = paste0(outdir, "/", argv$outprefix, "_")

# Load in files
seqtab.nochim <- readRDS(RDS_path)
otutab <- t(seqtab.nochim)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, paste0(refdir,refDB), multithread=TRUE, tryRC=TRUE)
taxa <- addSpecies(taxa, paste0(refDB, spref))

write.table(taxa, file = paste0(outprefix,"/taxa.tsv"), quote=FALSE)
saveRDS(taxa, paste0(outprefix,"/taxa.Rds"))

# Create abundance table with tax assignments and counts - merged by ASVs (rownames)
abundance <- merge(taxa, otutab, by = "row.names", all - TRUE)
write.csv(abundance, file = paste0(outprefix, "abundance.csv"), quote = FALSE)

# The following line replaces any instances of NNNNNNNNNN with nothing
# Use this if paired reads didn't merge in dada2, and you used justConcatenate = TRUE
# abundance$Row.names <- gsub("NNNNNNNNNN", "", abundance$Row.names)
