# This folder is for de-multiplexing reads by amplicon, removing primer sequences, and trimming for quality

## Demulitplexing
Run cutadapt in 3 steps. First, separate reverse primers HCO2198 and SSUR. Then, separate HCO2198 reads into forward primers LCO1490 and mICOintF. Then, go back and check unmatched read files for missed reads that the reverse primers missed.