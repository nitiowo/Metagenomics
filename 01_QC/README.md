# Quality control

Run FastQC on all raw reads, then summarize with MultiQC.
`sampsheet.sh` (in `00_prep/`) generates the sample list first.

- `run_qc.job` — SGE job that runs FastQC then MultiQC in sequence

## Configuration

Rename `QC_scripts/config.sh.example` to `config.sh` and edit parameters before running. `config.sh` is gitignored.