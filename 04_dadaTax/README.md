# Taxonomic assignment (DADA2)

Used `assignTaxonomy` + `addSpecies` against custom reference databases.
See `link to DB curation repo` for information on custom database creation.
Reference databases in `data/` directory.

- `dada_tax_assign.R` — takes a seqtab RDS and ref fastas, outputs taxa table
- `dada_tax.job` — SGE submission script