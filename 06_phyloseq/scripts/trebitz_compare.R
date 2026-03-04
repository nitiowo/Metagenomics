# trebitz_compare.R
# Compare detected taxa against Trebitz Great Lakes lists

source("setup.R")

trebitz <- load_trebitz()

outdir <- file.path(output_root, "trebitz_compare")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL
compare_rank <- "Genus"

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)
marker_names <- names(ps_filt)

# ---- Per-Lake Unexpected Taxa ----
per_lake <- imap_dfr(ps_filt, function(ps, mk) {
  compare_to_trebitz(ps, trebitz, rank = compare_rank) %>%
    mutate(Marker = mk)
})

save_stats(per_lake,
           file.path(outdir, "stats", "trebitz_per_lake_comparison"),
           caption = "Per-lake eDNA vs Trebitz comparison")

# ---- Overall Unexpected Taxa ----
overall <- per_lake %>%
  group_by(Marker) %>%
  summarise(n_edna_total = sum(n_edna),
            n_trebitz_total = sum(n_trebitz),
            n_shared_total = sum(n_shared),
            n_edna_only_total = sum(n_edna_only),
            .groups = "drop")

save_stats(overall,
           file.path(outdir, "stats", "trebitz_overall_summary"),
           caption = "Overall eDNA vs Trebitz overlap by marker")

# ---- Richness Comparison ----
richness_comp <- per_lake %>%
  select(Marker, Lake, n_edna, n_trebitz) %>%
  pivot_longer(cols = c(n_edna, n_trebitz),
               names_to = "Source", values_to = "Richness") %>%
  mutate(Source = recode(Source,
                         "n_edna" = "eDNA", "n_trebitz" = "Trebitz"))

p_rich <- ggplot(richness_comp,
                 aes(x = Lake, y = Richness, fill = Source)) +
  geom_col(position = "dodge", alpha = 0.8) +
  facet_wrap(~ Marker, scales = "free_y") +
  scale_fill_manual(values = c("eDNA" = "steelblue",
                                "Trebitz" = "coral")) +
  theme_minimal(base_size = 10) +
  labs(title = paste("Richness comparison:", compare_rank),
       y = paste(compare_rank, "richness"))
save_plot(p_rich, file.path(outdir, "figures",
                             "trebitz_richness_comparison.pdf"),
          width = 14, height = 8)

# ---- Scatterplot markers vs Trebitz Richness ----
scatter_df <- per_lake %>%
  select(Marker, Lake, n_edna, n_trebitz)

p_scat <- ggplot(scatter_df,
                  aes(x = n_trebitz, y = n_edna, color = Lake)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  facet_wrap(~ Marker) +
  scale_color_manual(values = lake_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "eDNA vs Trebitz richness per lake",
       x = "Trebitz genera", y = "eDNA genera")
save_plot(p_scat, file.path(outdir, "figures",
                             "trebitz_richness_scatter.pdf"),
          width = 12, height = 8)
