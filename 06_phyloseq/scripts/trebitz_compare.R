# trebitz_compare.R
# Compare detected taxa against Trebitz Great Lakes lists

source("setup.R")

outdir <- file.path(output_root, "trebitz_compare")

# ---- Control ----
use_mol_list  <- ps_markers           # markers only
compare_ranks <- c("Family")          # default Family only
use_lakes     <- NULL                 # NULL = all lakes

# trebitz_file is set in config.R

lake_cols <- c(Superior = "Superior", Michigan = "Michigan",
               Huron = "Huron", Erie = "Erie", Ontario = "Ontario")

# Unique non-NA taxon names from a ps object
get_det_taxa <- function(ps, rank, lake = NULL) {
  ps_agg <- agg_rank(ps, rank) %>% to_pa()
  if (!is.null(lake)) {
    sd  <- data.frame(sample_data(ps_agg))
    keep <- rownames(sd)[sd$Lake == lake]
    if (length(keep) == 0) return(character(0))
    ps_agg <- prune_samples(keep, ps_agg)
    ps_agg <- prune_taxa(taxa_sums(ps_agg) > 0, ps_agg)
  }
  tt <- data.frame(tax_table(ps_agg), stringsAsFactors = FALSE)
  unique(na.omit(tt[[rank]]))
}

# Unique non-NA taxon names from Trebitz list
get_treb_taxa <- function(treb_raw, rank, lake = NULL) {
  if (!(rank %in% colnames(treb_raw))) return(character(0))
  if (is.null(lake)) {
    return(unique(na.omit(treb_raw[[rank]])))
  }
  col <- lake_cols[[lake]]
  if (!(col %in% colnames(treb_raw))) return(character(0))
  treb_raw %>%
    filter(.data[[col]] == "X") %>%
    pull(.data[[rank]]) %>%
    na.omit() %>%
    unique()
}

# ---- Load Trebitz ----
raw <- NULL
if (file.exists(trebitz_file)) {
  raw <- read.csv(trebitz_file, stringsAsFactors = FALSE, check.names = FALSE)
  raw <- raw[, colnames(raw) != "" & !is.na(colnames(raw)), drop = FALSE]
  for (lk in names(lake_cols)) {
    col <- lake_cols[[lk]]
    if (col %in% colnames(raw)) {
      n <- sum(raw[[col]] == "X", na.rm = TRUE)
      cat("  ", lk, ":", n, "records\n")
    }
  }
} else {
  cat("Trebitz file not found:", trebitz_file, "\n")
}

active_lakes <- if (!is.null(use_lakes)) intersect(use_lakes, lake_order) else lake_order
# ---- Detection Rates vs Trebitz List ----
# Calculate % of Trebitz taxa detected by each method and combined molecular

scope_list <- c(setNames(as.list(active_lakes), active_lakes), list(Overall = NULL))
method_list <- c(use_mol_list, list(Morphology = ps_morph))

sens_df <- bind_rows(lapply(names(scope_list), function(lk_nm) {
  lk <- scope_list[[lk_nm]]
  rows <- bind_rows(lapply(names(method_list), function(m) {
    ps <- method_list[[m]]
    bind_rows(lapply(compare_ranks, function(r) {
      treb <- get_treb_taxa(raw, r, lk)
      det  <- get_det_taxa(ps, r, lk)
      n_ov <- length(intersect(det, treb))
      tibble(
        Lake          = lk_nm,
        Marker        = m,
        Rank          = r,
        n_trebitz     = length(treb),
        n_detected    = length(det),
        n_shared      = n_ov,
        n_novel       = length(setdiff(det, treb)),
        detection_pct = if (length(treb) > 0) round(100 * n_ov / length(treb), 1) else NA_real_
      )
    }))
  }))
  # Combined molecular
  comb <- bind_rows(lapply(compare_ranks, function(r) {
    treb <- get_treb_taxa(raw, r, lk)
    det  <- unique(unlist(lapply(use_mol_list, get_det_taxa, rank = r, lake = lk)))
    n_ov <- length(intersect(det, treb))
    tibble(
      Lake          = lk_nm,
      Marker        = "Combined_Mol",
      Rank          = r,
      n_trebitz     = length(treb),
      n_detected    = length(det),
      n_shared      = n_ov,
      n_novel       = length(setdiff(det, treb)),
      detection_pct = if (length(treb) > 0) round(100 * n_ov / length(treb), 1) else NA_real_
    )
  }))
  bind_rows(rows, comb)
}))

sens_df$Rank   <- factor(sens_df$Rank,   levels = compare_ranks, ordered = TRUE)
sens_df$Marker <- factor(sens_df$Marker, levels = c(names(use_mol_list), "Morphology", "Combined_Mol"))
sens_df$Lake   <- factor(sens_df$Lake,   levels = c(active_lakes, "Overall"))

# save the stats first
save_stats(sens_df,
           file.path(outdir, "stats", "detection_by_rank_lake"),
           caption = "Detection rates vs Trebitz published list (per rank, lake, and method)")

method_pal <- c(marker_colors, Combined_Mol = "#333333")

# Plot
p_sens <- ggplot(
    filter(sens_df, !is.na(detection_pct)),
    aes(x = Lake, y = detection_pct, fill = Marker)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Rank, ncol = 1) +
  scale_fill_manual(values = method_pal, name = "Method") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = plot_base_size) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text   = element_text(face = "bold")) +
  labs(title = "Detection rate vs Trebitz published list",
       x = NULL, y = "% of Trebitz taxa detected")

# Save plot
save_plot(p_sens,
          file.path(outdir, "figures", "detection_by_rank_lake.pdf"),
          width = 12, height = 10)




# ---- Three-way Comparison: Trebitz vs Morphology vs Combined Molecular ----

cat_levels <- c("treb_only", "treb_morph", "treb_mol", "all_three",
                "morph_only", "morph_mol",  "mol_only")
cat_labels <- c("Trebitz only", "Trebitz + Morph", "Trebitz + Mol",
                "All three", "Morph only", "Morph + Mol", "Mol only")
cat_colors <- c(
  "Trebitz only"    = "#d9d9d9",
  "Trebitz + Morph" = "#a6cee3",
  "Trebitz + Mol"   = "#fb9a99",
  "All three"       = "#b2df8a",
  "Morph only"      = "#1f78b4",
  "Morph + Mol"     = "#6a51a3",
  "Mol only"        = "#e31a1c"
)

# Check that this is correct - code runs but may not do
threeway_df <- bind_rows(lapply(names(scope_list), function(lk_nm) {
  lk <- scope_list[[lk_nm]]
  bind_rows(lapply(compare_ranks, function(r) {
    treb  <- get_treb_taxa(raw, r, lk)
    morph <- get_det_taxa(ps_morph, r, lk)
    mol   <- unique(unlist(lapply(use_mol_list, get_det_taxa, rank = r, lake = lk)))
    all_taxa <- unique(c(treb, morph, mol))
    in_treb  <- all_taxa %in% treb
    in_morph <- all_taxa %in% morph
    in_mol   <- all_taxa %in% mol
    tibble(
      Lake       = lk_nm, Rank = r,
      treb_only  = sum( in_treb & !in_morph & !in_mol),
      morph_only = sum(!in_treb &  in_morph & !in_mol),
      mol_only   = sum(!in_treb & !in_morph &  in_mol),
      treb_morph = sum( in_treb &  in_morph & !in_mol),
      treb_mol   = sum( in_treb & !in_morph &  in_mol),
      morph_mol  = sum(!in_treb &  in_morph &  in_mol),
      all_three  = sum( in_treb &  in_morph &  in_mol),
      n_treb = length(treb), n_morph = length(morph), n_mol = length(mol)
    )
  }))
}))

threeway_df$Rank <- factor(threeway_df$Rank, levels = compare_ranks, ordered = TRUE)
threeway_df$Lake <- factor(threeway_df$Lake, levels = c(active_lakes, "Overall"))

save_stats(threeway_df,
           file.path(outdir, "stats", "threeway_overlap"),
           caption = "Three-way overlap: Trebitz / Morphology / Combined Molecular")

tw_long <- threeway_df %>%
  pivot_longer(cols = all_of(cat_levels),
               names_to = "Category", values_to = "n_taxa") %>%
  mutate(Category = factor(
    recode(Category,
      treb_only  = "Trebitz only",
      morph_only = "Morph only",
      mol_only   = "Mol only",
      treb_morph = "Trebitz + Morph",
      treb_mol   = "Trebitz + Mol",
      morph_mol  = "Morph + Mol",
      all_three  = "All three"
    ), levels = cat_labels))

p_three <- ggplot(tw_long, aes(x = Lake, y = n_taxa, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", colour = "white", linewidth = 0.2) +
  facet_wrap(~ Rank, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cat_colors, name = "Detection category") +
  theme_minimal(base_size = plot_base_size) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text   = element_text(face = "bold")) +
  labs(title = "Three-way comparison: Trebitz / Morphology / Molecular",
       x = NULL, y = "Number of taxa")

save_plot(p_three,
          file.path(outdir, "figures", "threeway_overlap_by_lake.pdf"),
          width = 12, height = 10)

