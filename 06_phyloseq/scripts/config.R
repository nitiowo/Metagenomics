# config.R
# Master configuration — edit this file to run the pipeline from any location
# All other scripts read these variables via setup.R (which sources this file)

# ---- Pipeline root ----
# Directory containing the analysis scripts (this file lives here)
pipeline_dir <- "/Volumes/Samsung_1TB/Zooplankton/git_DBs/Github_config/testing/06_phyloseq"

# ---- Output root ----
# All output subdirectories (figures/, stats/) are created under here
output_root <- "output"

# ---- Input RDS paths ----
rds_folmer <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/setup_output/folmer_ps.RDS"
rds_leray  <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/setup_output/leray_ps.RDS"
rds_18S    <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/setup_output/ssu_ps.RDS"
rds_morph  <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/setup_output/morph_ps.RDS"

# ---- Published reference list (Trebitz) ----
trebitz_file <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/data/trebitz_lists/Trebitz_Zoops_2026_overall_plus_benthos_taxfixed.csv"
