# geographic.R
# Geographic patterns: map, IBD, distance-decay, centroids, dbMEMs, NMDS map

source("setup.R")

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

# ---- Distance-Decay Curves ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  sdf_mk <- filter_geo_metadata(ps)
  ps_sub <- prune_samples(rownames(sdf_mk), ps)
  
  comm_dist <- phyloseq::distance(ps_sub, method = "jaccard")
  coords <- sdf_mk %>% select(Longitude, Latitude) %>% as.matrix()
  geo_m <- geosphere::distm(coords, fun = geosphere::distHaversine)
  
  dd_df <- tibble(
    geo_km = as.vector(as.dist(geo_m)) / 1000,
    comm_dist = as.vector(comm_dist)
  )
  
  p_dd <- ggplot(dd_df, aes(x = geo_km, y = comm_dist)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", color = "red") +
    theme_minimal(base_size = 10) +
    labs(title = paste("Distance-decay:", mk),
         x = "Geographic distance (km)",
         y = "Community dissimilarity (Jaccard)")
  save_plot(p_dd, file.path(outdir, "figures",
                             paste0("distance_decay_", mk, ".pdf")))
}

# ---- Centroid Distances Between Lakes ----
centroids <- sdf %>%
  group_by(Lake) %>%
  summarise(Lon = mean(Longitude), Lat = mean(Latitude), .groups = "drop")

centroid_pairs <- combn(centroids$Lake, 2, simplify = FALSE)
cent_dists <- map_dfr(centroid_pairs, function(pr) {
  c1 <- centroids %>% filter(Lake == pr[1])
  c2 <- centroids %>% filter(Lake == pr[2])
  d <- geosphere::distHaversine(c(c1$Lon, c1$Lat), c(c2$Lon, c2$Lat))
  tibble(Lake_A = pr[1], Lake_B = pr[2], dist_km = d / 1000)
})

save_stats(cent_dists,
           file.path(outdir, "stats", "centroid_distances_between_lakes"),
           caption = "Geographic distance between lake centroids (km)")

# ---- dbMEMs Per Marker ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  sdf_mk <- filter_geo_metadata(ps)
  ps_sub <- prune_samples(rownames(sdf_mk), ps)
  
  coords <- sdf_mk %>% select(Longitude, Latitude) %>% as.matrix()
  geo_d <- geosphere::distm(coords, fun = geosphere::distHaversine)
  
  dbmem_res <- tryCatch(
    adespatial::dbmem(as.dist(geo_d), silent = TRUE),
    error = function(e) NULL
  )
  if (is.null(dbmem_res)) next
  
  # RDA with positive dbMEMs
  otu <- as.data.frame(t(otu_table(ps_sub)))
  if (taxa_are_rows(ps_sub)) otu <- as.data.frame(otu_table(ps_sub)) %>% t() %>% as.data.frame()
  
  hel <- vegan::decostand(otu, method = "hellinger")
  mem_df <- as.data.frame(dbmem_res)
  
  rda_res <- tryCatch(
    vegan::rda(hel ~ ., data = mem_df),
    error = function(e) NULL
  )
  if (is.null(rda_res)) next
  
  rda_test <- vegan::anova.cca(rda_res, permutations = 999)
  save_stats(tibble(Marker = mk,
                    Variance_explained = RsquareAdj(rda_res)$adj.r.squared,
                    F_stat = rda_test$F[1],
                    p_value = rda_test$`Pr(>F)`[1]),
             file.path(outdir, "stats",
                        paste0("dbmem_rda_", mk)),
             caption = paste("dbMEM RDA:", mk))
}

# ---- NMDS Map Overlay ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  sdf_mk <- filter_geo_metadata(ps)
  ps_sub <- prune_samples(rownames(sdf_mk), ps)
  
  dist_mat <- phyloseq::distance(ps_sub, method = "jaccard")
  ord <- ordinate(ps_sub, "NMDS", distance = dist_mat, trymax = 50)
  
  ord_df <- data.frame(
    NMDS1 = ord$points[, 1],
    NMDS2 = ord$points[, 2],
    Longitude = sdf_mk$Longitude,
    Latitude = sdf_mk$Latitude,
    Lake = sdf_mk$Lake
  )
  
  p_nm <- ggplot(ord_df, aes(x = Longitude, y = Latitude,
                              color = NMDS1, shape = Lake)) +
    geom_point(size = 3.5) +
    scale_color_viridis_c() +
    theme_minimal(base_size = 10) +
    labs(title = paste("NMDS axis 1 mapped spatially:", mk))
  save_plot(p_nm, file.path(outdir, "figures",
                             paste0("nmds_map_overlay_", mk, ".pdf")),
            width = 10, height = 7)
}
