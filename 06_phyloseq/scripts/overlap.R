# overlap.R
# Taxa overlap: Venn diagrams, pairwise matrix, UpSet plots

source("setup.R")

outdir <- file.path(output_root, "overlap")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL
overlap_rank <- "Genus"

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)
marker_names <- names(ps_filt)

# ---- Per-Marker Venn Diagrams ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  taxa_sets <- get_taxa_set(ps, rank = overlap_rank)
  
  if (length(taxa_sets) < 2 || length(taxa_sets) > 5) next
  
  cols <- lake_colors[names(taxa_sets)]
  
  venn.plot <- VennDiagram::venn.diagram(
    x = taxa_sets,
    fill = cols,
    alpha = 0.5,
    filename = NULL,
    main = paste(mk, "-", overlap_rank, "overlap by lake")
  )
  
  pdf(file.path(outdir, "figures",
                paste0("venn_", mk, "_", overlap_rank, ".pdf")),
      width = 8, height = 8)
  grid::grid.draw(venn.plot)
  dev.off()
}

# ---- Pairwise Overlap Matrix ----
pw_overlap <- imap_dfr(ps_filt, function(ps, mk) {
  taxa_sets <- get_taxa_set(ps, rank = overlap_rank)
  lakes <- names(taxa_sets)
  pairs <- combn(lakes, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pr) {
    shared <- length(intersect(taxa_sets[[pr[1]]], taxa_sets[[pr[2]]]))
    union_n <- length(union(taxa_sets[[pr[1]]], taxa_sets[[pr[2]]]))
    tibble(Marker = mk, Lake_A = pr[1], Lake_B = pr[2],
           Shared = shared, Union = union_n,
           Jaccard = shared / union_n)
  })
})

save_stats(pw_overlap,
           file.path(outdir, "stats",
                      paste0("pairwise_overlap_", overlap_rank)),
           caption = paste("Pairwise taxa overlap:", overlap_rank))

# ---- Between-Marker UpSet ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]
  upset_mat <- make_upset_matrix(ps, rank = overlap_rank)
  
  pdf(file.path(outdir, "figures",
                paste0("upset_", mk, "_", overlap_rank, ".pdf")),
      width = 10, height = 6)
  print(UpSetR::upset(upset_mat, order.by = "freq",
                      main.bar.color = "steelblue",
                      sets.bar.color = "grey40",
                      text.scale = 1.2))
  dev.off()
}

# Between-marker overlap: taxa shared across markers at each lake
cross_marker_sets <- map(unique(unlist(map(ps_filt, ~ {
  data.frame(sample_data(.x))$Lake
}))), function(lk) {
  imap(ps_filt, function(ps, mk) {
    sdf <- data.frame(sample_data(ps))
    samps <- rownames(sdf[sdf$Lake == lk, ])
    if (length(samps) == 0) return(character(0))
    ps_sub <- prune_samples(samps, ps)
    ps_agg <- agg_rank(ps_sub, overlap_rank)
    taxa_names(ps_agg)[taxa_sums(ps_agg) > 0]
  })
}) %>% setNames(unique(unlist(map(ps_filt, ~ {
  data.frame(sample_data(.x))$Lake
}))))

# Cross-marker UpSet per lake
for (lk in names(cross_marker_sets)) {
  sets <- cross_marker_sets[[lk]]
  if (length(sets) < 2) next
  all_taxa <- unique(unlist(sets))
  mat <- data.frame(row.names = all_taxa)
  for (mk in names(sets)) {
    mat[[mk]] <- as.integer(all_taxa %in% sets[[mk]])
  }
  
  pdf(file.path(outdir, "figures",
                paste0("upset_cross_marker_", lk, ".pdf")),
      width = 10, height = 6)
  print(UpSetR::upset(mat, order.by = "freq",
                      main.bar.color = "steelblue",
                      text.scale = 1.2))
  dev.off()
}
