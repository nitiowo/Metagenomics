# ASV inference (DADA2)

Run DADA2 on each primer set separately.
- `dada2_asv.R` — paired-end pipeline
- `dada2_asv_SE.R` — single-end variant (used for LCO1490 which didn't merge well)
- `dada_asv.job` — SGE submission script

## Configuration

Rename `dada_asv_scripts/config.sh.example` to `config.sh` and edit parameters before running. `config.sh` is gitignored.