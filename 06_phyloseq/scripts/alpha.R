# alpha.R
# Alpha diversity: compute and plot metrics

source("setup.R")

outdir <- file.path(output_root, "alpha")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL
use_metrics <- c("Observed", "InvSimpson")

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)

# ---- Compute a div across all markers----
alpha_all <- compute_alpha_all(ps_filt, use_metrics, lake_order)

alpha_long <- alpha_all %>%
  pivot_longer(cols = all_of(use_metrics),
               names_to = "Metric", values_to = "Value")

# ---- Boxplot: Between Markers ----
p_mk <- ggplot(alpha_long, aes(x = Marker, y = Value, fill = Marker)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = marker_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Alpha diversity - between markers",
       y = "Diversity value") +
  theme(legend.position = "bottom")
save_plot(p_mk, file.path(outdir, "figures",
                           "alpha_boxplot_between_markers.pdf"))

# ---- Boxplot: Lake by Marker ----
p_lk <- ggplot(alpha_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_grid(Metric ~ Marker, scales = "free_y") +
  scale_fill_manual(values = lake_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Alpha diversity - between lakes (per marker)",
       y = "Diversity value") +
  theme(legend.position = "bottom")
save_plot(p_lk, file.path(outdir, "figures",
                           "alpha_boxplot_lake_by_marker.pdf"),
          width = 16, height = 10)
