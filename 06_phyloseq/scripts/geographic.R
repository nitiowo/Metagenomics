outdir <- file.path(output_root, "geographic")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)
marker_names <- names(ps_filt)

# ---- Station Map ----
sdf <- filter_geo_metadata(ps_filt[[1]])

p_map <- ggplot(sdf, aes(x = Longitude, y = Latitude, color = Lake)) +
  geom_point(size = 3) +
  scale_color_manual(values = lake_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Sampling stations")
save_plot(p_map, file.path(outdir, "figures", "station_map.pdf"),
          width = 10, height = 7)

# ---- Mantel IBD Per Marker ----
ibd_results <- map_dfr(marker_names, function(mk) {
  run_ibd(ps_filt[[mk]], dist_method = "jaccard") %>%
    mutate(Marker = mk)
})

save_stats(ibd_results,
           file.path(outdir, "stats", "mantel_IBD_per_marker"),
           caption = "Mantel test: isolation by distance (Jaccard)")
