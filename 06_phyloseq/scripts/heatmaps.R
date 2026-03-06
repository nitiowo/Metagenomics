# heatmaps.R
# Top-N heatmaps per marker with taxa subsetting

source("setup.R")

outdir <- file.path(output_root, "heatmaps")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL
top_n <- 30
heatmap_rank <- "Genus"

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)
marker_names <- names(ps_filt)

# ---- Top-N Heatmaps ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  plot_top_heatmap(ps, marker_name = mk, top_n = top_n,
                   rank = heatmap_rank,
                   outfile = file.path(outdir, "figures",
                                        paste0("heatmap_top", top_n, "_",
                                               mk, ".pdf")))
}

# ---- Subsetted Heatmaps (Phylum-Level) ----
focal_phyla <- c("Arthropoda", "Chordata", "Rotifera")

for (phylum in focal_phyla) {
  for (mk in marker_names) {
    ps <- ps_filt[[mk]]
    ps_sub <- subset_focal_taxon(ps, "Phylum", phylum)
    if (is.null(ps_sub) || ntaxa(ps_sub) == 0) next
    
    plot_top_heatmap(ps_sub, marker_name = paste(mk, phylum, sep = "_"),
                     top_n = min(top_n, ntaxa(agg_rank(ps_sub, heatmap_rank))),
                     rank = heatmap_rank,
                     outfile = file.path(outdir, "figures",
                                          paste0("heatmap_", phylum, "_",
                                                 mk, ".pdf")))
  }
}
