# Great Lakes Zooplankton Metabarcoding

Scripts for processing zooplankton metabarcoding samples from the Great Lakes.
Three primer sets were used (Folmer COI, Leray COI, 18S rRNA) alongside
morphological identification.

## Pipeline steps

| Dir | Step |
|-----|------|
| `00_prep` | Generate sample lists from raw fastq files |
| `01_QC` | FastQC + MultiQC on raw reads |
| `02_cutadapt` | Demultiplex by primer (and primer combos) and trim adapters, plus post-trim QC |
| `03_dadaASV` | DADA2 ASV inference |
| `04_dadaTax` | DADA2 taxonomic assignment |
| `05_taxize` | Taxonomic cleanup with taxize |
| `06_phyloseq` | Community analysis in phyloseq |
| `data` | Raw data, metadata, and reference database files |

Raw data and large outputs are gitignored. See `.gitignore` for details.