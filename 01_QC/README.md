# Quality control

Run FastQC on all raw reads, then summarize with MultiQC.
`sampsheet.sh` (in `00_prep/`) generates the sample list first.

- `fastQC.sh` — runs fastqc on everything in the raw data dir
- `multiqc.sh` — aggregates fastqc output
- `fullQC.job` — SGE job that runs both