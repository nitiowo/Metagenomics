# differential.R
# Differential analysis: SIMPER + indicator species

source("setup.R")

outdir <- file.path(output_root, "differential")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)
marker_names <- names(ps_filt)

# ---- SIMPER Per Marker ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  sim <- run_simper_analysis(ps, group_var = "Lake")
  if (is.null(sim)) next

  save_stats(sim, file.path(outdir, "stats",
                             paste0("simper_top_contributors_", mk)),
             caption = paste("SIMPER: top contributors", mk))
}

# ---- Indicator Species Per Marker ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  ind <- run_indicator(ps, group_var = "Lake")

  if (nrow(ind) == 0) next
  save_stats(ind, file.path(outdir, "stats",
                             paste0("indicator_species_", mk)),
             caption = paste("Indicator species:", mk))
}
