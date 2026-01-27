# This is inspired by https://benjjneb.github.io/dada2/tutorial_1_8.html
# Usage: dada2_asv.R --input </inputdir/path> --output </outputdir/path> --outprefix "prefix"

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.11")
# install.packages('argparser')


library(dada2); packageVersion("dada2")
library(argparser, quietly = TRUE)

p <- arg_parser("DADA2 in R")

p = add_argument(p, "--input", help = "Enter input folder", type = 'character')
p = add_argument(p, "--output", help = "Enter output folder", type = 'character')
p = add_argument(p, "--outprefix", help = "Enter output prefix for filenames", type = 'character')

# Parse the command line arguments
argv <- parse_args(p)

# SCRIPT MODE
path = argv$input
outdir = argv$output
outprefix = paste0(outdir, "/", argv$outprefix, "_")

# INTERACTIVE MODE
# path = '/Volumes/Seagate_2/Bioinformatics/Zooplankton/zoop_backup/subset/cutadapt/cutadapt_demux/mICOint'
# outdir = '/Volumes/Seagate_2/Bioinformatics/Zooplankton/zoop_backup/subset/cutadapt/cutadapt_demux/mICOint/dada_output/'
# outprefix = "mICOint"

files <- list.files(path)
fastqs <- files[grepl(".fastq.gz$", files)] # gz
# fastq_files <- files[grepl(".fastq", files)]


######################################################################################################
# Setting naming convention for F/R input fastq files
######################################################################################################
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), 
                       function(x) paste0("Z", x[1]))  # This part adds "Z" in front of each sample number
# Change based on sample names

# This splits by "_" and takes the first field as sample names
# sapply(strsplit(basename(fnFs), "_"), `[`, 1)


####################################################################
# Visualizing the quality profiles of the forward and reverse reads
####################################################################
pdf(paste0(outprefix, "dada2_plotQualityProfile.pdf"), onefile = TRUE)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
dev.off()


#############################################################
# Perform filtering and trimming on forward and reverse reads
##############################################################
# Place under filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Modify and add new parameter based on your data.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(272,218), maxN = 0, maxEE = c(2,3), truncQ = 5, 
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE, minLen = 40)


# Examine quality profiles of filtered reads
pdf(paste0(outprefix, "QualityProfile.filt_plot.pdf"), onefile = T)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
dev.off()


#####################################################
# Learn the Error Rates for forward and reverse reads
#####################################################

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

# Plot estimated error as sanity check
pdf(paste0(outprefix, "ErrorsRates_F.pdf"), onefile = TRUE)
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(paste0(outprefix, "ErrorsRates_R.pdf"), onefile = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()


################
# Dereplication
################
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


###################
# Sample Inference
###################
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


####################
# Merge paired reads
####################
# Make sure your amplicon length allows 20+ bp of overlap between F and R reads
# If not, try justConcatenate = TRUE (This just joins the two reads with 10 Ns in between). Code to remove the Ns is provided below.
# Or just use forward reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
# head(mergers[[1]])


##########################
# Construct sequence table
###########################

seqtab <- makeSequenceTable(mergers)
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
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, paste0(outprefix, "track_reads.txt"), sep = "\t", quote = FALSE)



