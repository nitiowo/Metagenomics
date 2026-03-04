# varpart.R
# Variance partitoning using spatial data vs Lake

source("setup.R")

outdir <- file.path(output_root, "varpart")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL

# ---- Per-Marker Varpart ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)
marker_names <- names(ps_filt)

vp_results <- list()

for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  vp <- run_varpart(ps)
  if (is.null(vp)) next
  vp_results[[mk]] <- vp

  pdf(file.path(outdir, "figures", paste0("varpart_venn_", mk, ".pdf")),
      width = 6, height = 6)
  plot(vp, bg = c("steelblue", "tomato"), Xnames = c("Lake", "Station"),
       main = mk)
  dev.off()
}

# Fractions table
frac_table <- imap_dfr(vp_results, function(vp, mk) {
  fr <- vp$part$indfract
  tibble(Marker = mk,
         Lake_only = fr$Adj.R.squared[1],
         Station_only = fr$Adj.R.squared[2],
         Shared = fr$Adj.R.squared[3],
         Residual = fr$Adj.R.squared[4])
})

save_stats(frac_table,
           file.path(outdir, "stats", "varpart_fractions_per_marker"),
           caption = "Variance partitioning: adjusted R-squared fractions")

# ---- Combined Markers Varpart ----
for (ps_comb_name in c("ps_markers_combined", "ps_coi_combined")) {
  ps_c <- get(ps_comb_name, envir = .GlobalEnv)
  label <- gsub("ps_|_combined", "", ps_comb_name)

  vp_c <- run_varpart(ps_c)
  if (is.null(vp_c)) next

  pdf(file.path(outdir, "figures",
                paste0("varpart_venn_combined_", label, ".pdf")),
      width = 6, height = 6)
  plot(vp_c, bg = c("steelblue", "tomato"),
       Xnames = c("Lake", "Station"),
       main = paste("Combined:", label))
  dev.off()

  fr_c <- vp_c$part$indfract
  frac_c <- tibble(Dataset = label,
                   Lake_only = fr_c$Adj.R.squared[1],
                   Station_only = fr_c$Adj.R.squared[2],
                   Shared = fr_c$Adj.R.squared[3],
                   Residual = fr_c$Adj.R.squared[4])
  save_stats(frac_c,
             file.path(outdir, "stats",
                        paste0("varpart_fractions_combined_", label)),
             caption = paste("Varpart fractions: combined", label))
}
