# This is inspired by https://benjjneb.github.io/dada2/tutorial_1_8.html
# Usage: dada2_asv.R --input </inputdir/path> --output </outputdir/path> --outprefix "prefix"

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.11")
# install.packages('argparser')

# This scripts only uses forward reads and skips merging before ASV generation

library(dada2); packageVersion("dada2")
library(argparser, quietly = TRUE)

p <- arg_parser("DADA2 in R")

p = add_argument(p, "--input", help = "Enter input folder", type = 'character')
p = add_argument(p, "--output", help = "Enter output folder", type = 'character')
p = add_argument(p, "--outprefix", help = "Enter output prefix for filenames", type = 'character')
p = add_argument(p, "--trunclen", help = "truncLen for forward reads (single value)", type = 'integer', default = 273)
p = add_argument(p, "--maxee", help = "maxEE for forward reads (single value)", type = 'numeric', default = 3)
p = add_argument(p, "--maxn", help = "Max Ns allowed", type = 'integer', default = 0)
p = add_argument(p, "--minlen", help = "Minimum read length after filtering", type = 'integer', default = 40)

# Parse the command line arguments
argv <- parse_args(p)

# SCRIPT MODE
path = argv$input
outdir = argv$output
outprefix = paste0(outdir, "/", argv$outprefix, "_")

# Parse filtering params
truncLen_val <- argv$trunclen
maxEE_val <- argv$maxee
maxN_val <- argv$maxn
minLen_val <- argv$minlen

files <- list.files(path)
fastqs <- files[grepl(".fastq.gz$", files)] # gz
# fastq_files <- files[grepl(".fastq", files)]

dir.create(outdir)

######################################################
# Setting naming convention for F/R input fastq files
######################################################
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
# fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), 
                       function(x) paste0("Z", x[1]))  # This part adds "Z" in front of each sample number
# Change based on sample names

# This splits by "_" and takes the first field as sample names
# sapply(strsplit(basename(fnFs), "_"), `[`, 1)


####################################################################
# Visualizing the quality profiles of the forward and reverse reads
####################################################################
pdf(paste0(outprefix, "dada2_plotQualityProfile_fwd.pdf"), onefile = TRUE)
plotQualityProfile(fnFs[1:2])
dev.off()

# pdf(paste0(outprefix, "dada2_plotQualityProfile_rev.pdf"), onefile = TRUE)
# plotQualityProfile(fnRs[1:2])
# dev.off()

#############################################################
# Perform filtering and trimming on forward and reverse reads
##############################################################
# Place under filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
# names(filtRs) <- sample.names

# Choose trunclen based on quality profiles
out <- filterAndTrim(fnFs, filtFs, truncLen = truncLen_val, maxN = maxN_val, maxEE = maxEE_val, 
                     #truncQ=5, truncQ not recommended
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE, minLen = minLen_val)


# Examine quality profiles of filtered reads
pdf(paste0(outprefix, "QualityProfile.filt_plot_fwd.pdf"), onefile = T)
plotQualityProfile(filtFs[1:2])
dev.off()

# pdf(paste0(outprefix, "QualityProfile.filt_plot_rev.pdf"), onefile = T)
# plotQualityProfile(filtRs[1:2])
# dev.off()

#####################################################
# Learn the Error Rates for forward and reverse reads
#####################################################

errF <- learnErrors(filtFs, multithread = TRUE)

# errR <- learnErrors(filtRs, multithread = TRUE)

# Plot estimated error as sanity check
pdf(paste0(outprefix, "ErrorsRates_F.pdf"), onefile = TRUE)
plotErrors(errF, nominalQ = TRUE)
dev.off()

# pdf(paste0(outprefix, "ErrorsRates_R.pdf"), onefile = TRUE)
# plotErrors(errR, nominalQ = TRUE)
# dev.off()


################
# Dereplication
################
derepFs <- derepFastq(filtFs, verbose = TRUE)
# derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
# names(derepRs) <- sample.names


###################
# Sample Inference
###################
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
# dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


####################
# Merge paired reads
####################
# Make sure your amplicon length allows 20+ bp of overlap between F and R reads
# If not, try justConcatenate = TRUE (This just joins the two reads with 10 Ns in between). You can remove the Ns later
# But the recommended method is to just use forward reads

# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
# head(mergers[[1]])


##########################
# Construct sequence table
###########################

seqtab <- makeSequenceTable(dadaFs)
## Get dimensions
# dim(seqtab)

## Inspect distribution of sequence lengths
# table(nchar(getSequences(seqtab)))


##################
# Remove chimeras
##################
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)
# dim(seqtab.nochim)
# sum(seqtab.nochim)/sum(seqtab)


#################
#Save into Rds
#################
saveRDS(seqtab, paste0(outprefix, "seqtab.Rds"))
saveRDS(seqtab.nochim, paste0(outprefix, "seqtab.nochim.Rds"))


#Save seqtab as otu table
otutab <- t(seqtab.nochim)
write.table(otutab, file = paste0(outprefix, "otutab.tsv"), quote = FALSE)


#################################
# Track reads through the pipeline
#################################
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
track <- cbind(track, (track[,ncol(track)]/track[,1])*100)

colnames(track) <- c("input", "filtered", "denoisedF", "nonchim", "%_reads_remaining")
rownames(track) <- sample.names
head(track)
write.table(track, paste0(outprefix, "track_reads.txt"), sep = "\t",quote = FALSE)
