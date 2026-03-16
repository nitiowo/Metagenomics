# exploratory.R
# Dataset summaries: counts, unique taxa, % unassigned at each rank

source("setup.R")

outdir <- file.path(output_root, "exploratory")

# ---- Control ----
use_ps_list <- ps_markers   # molecular markers for tax assignment plots
top_n_taxa  <- 4            # top N named taxa colored per rank; rest = "Other"

# ---- Dataset Summary Tables ----

datasets <- c(ps_all_methods,
              list(`Combined markers` = ps_markers_combined,
                   `Combined all` = ps_all_combined))

summary_df <- imap_dfr(datasets, ~ summarise_ps(.x, .y, tax_ranks))
summary_df <- summary_df %>%
  select(dataset, samples, total_asvs, rank, unique_taxa,
         unassigned_asvs, pct_reads_unassigned)

save_stats(summary_df,
           file.path(outdir, "stats", "dataset_summary"),
           caption = "Dataset summary: taxa counts and % reads unassigned")

txt <- capture.output({
  for (nm in names(datasets))
    summarise_ps_print(datasets[[nm]], nm, tax_ranks)
})
cat(paste(txt, collapse = "\n"), "\n")
save_summary(txt, file.path(outdir, "stats", "dataset_summary.txt"))

# ---- Taxonomic Assignment Profile (% Abundance at Each Rank) ----

# Build per-marker x rank abundance fractions
abund_df <- imap_dfr(use_ps_list, function(ps, m) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  tt <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
  asv_abund <- rowSums(otu)

  imap_dfr(tax_ranks, function(rank, .i) {
    if (!(rank %in% colnames(tt))) return(tibble())
    # Species: use full binomial (Genus + epithet)
    if (rank == "Species" && "Genus" %in% colnames(tt)) {
      genus <- tt[["Genus"]]
      genus[is.na(genus) | genus == "" | genus == "NA"] <- NA
      sp    <- tt[["Species"]]
      sp[is.na(sp) | sp == "" | sp == "NA"] <- NA
      vals  <- ifelse(!is.na(genus) & !is.na(sp), paste(genus, sp), sp)
    } else {
      vals <- tt[[rank]]
      vals[is.na(vals) | vals == "" | vals == "NA"] <- NA
    }

    # Sum abundance for each taxon name
    agg <- tibble(Taxon = vals, abund = asv_abund) %>%
      mutate(Taxon = ifelse(is.na(Taxon), "Unassigned", Taxon)) %>%
      group_by(Taxon) %>%
      summarise(abund = sum(abund), .groups = "drop")

    total <- sum(agg$abund)
    agg %>% mutate(Marker = m, Rank = rank,
                   Pct = if (total > 0) abund / total else 0)
  })
})

abund_df$Rank   <- factor(abund_df$Rank, levels = tax_ranks, ordered = TRUE)
abund_df$Marker <- factor(abund_df$Marker, levels = names(use_ps_list))

# Per rank: keep top N taxa by total abundance, lump rest to "Other"
abund_plot <- abund_df %>%
  group_by(Rank) %>%
  group_modify(~ {
    named <- .x %>% filter(Taxon != "Unassigned")
    top <- named %>%
      group_by(Taxon) %>%
      summarise(tot = sum(abund), .groups = "drop") %>%
      arrange(desc(tot)) %>%
      slice_head(n = top_n_taxa) %>%
      pull(Taxon)
    .x %>% mutate(Taxon_plot = case_when(
      Taxon == "Unassigned" ~ "Unassigned",
      Taxon %in% top        ~ Taxon,
      TRUE                  ~ "Other"
    ))
  }) %>% ungroup() %>%
  group_by(Marker, Rank, Taxon_plot) %>%
  summarise(Pct = sum(Pct), .groups = "drop")

# per-rank top taxa get distinct colors
all_named <- abund_plot %>%
  filter(!Taxon_plot %in% c("Other", "Unassigned")) %>%
  distinct(Rank, Taxon_plot) %>%
  arrange(Rank)
n_named <- nrow(all_named)
named_cols <- colorRampPalette(brewer.pal(12, "Paired"))(max(n_named, 1))
names(named_cols) <- all_named$Taxon_plot

fill_pal <- c(named_cols, Other = "grey60", Unassigned = "grey85")

present_taxa <- unique(abund_plot$Taxon_plot)
fill_pal <- fill_pal[names(fill_pal) %in% present_taxa]

#Unassigned on top
taxa_stack <- setdiff(names(fill_pal), c("Other", "Unassigned"))
abund_plot$Taxon_plot <- factor(
  abund_plot$Taxon_plot,
  levels = rev(c(taxa_stack, "Other", "Unassigned")))

p_tax <- ggplot(abund_plot,
                aes(x = Marker, y = Pct * 100, fill = Taxon_plot)) +
  geom_bar(stat = "identity", position = "stack", colour = "white",
           linewidth = 0.2) +
  facet_grid(. ~ Rank, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fill_pal) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  theme_minimal(base_size = plot_base_size) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  ) +
  labs(
    title = "Taxonomic assignment profile (% abundance) by rank and marker",
    subtitle = paste("Top", top_n_taxa, "taxa colored per rank; remainder lumped to Other"),
    x = NULL, y = "% of total read abundance"
  )

# Assignment rate table
assign_rate <- abund_plot %>%
  filter(Taxon_plot != "Unassigned") %>%
  group_by(Marker, Rank) %>%
  summarise(pct_assigned = round(sum(Pct) * 100, 1), .groups = "drop") %>%
  arrange(Marker, Rank)

save_stats(assign_rate,
           file.path(outdir, "stats", "assignment_rates_abundance"),
           caption = "% read abundance assigned at each rank per marker")

# ---- ASV-count-based Assignment Rates ----

asv_pct_df <- imap_dfr(use_ps_list, function(ps, m) {
  tt <- data.frame(tax_table(ps), stringsAsFactors = FALSE)

  imap_dfr(tax_ranks, function(rank, .i) {
    if (!(rank %in% colnames(tt))) return(tibble())

    vals <- tt[[rank]]
    vals[is.na(vals) | vals == "" | vals == "NA"] <- NA

    counts <- table(vals, useNA = "always")
    names(counts)[is.na(names(counts))] <- "Unassigned"

    tibble(
      Marker = m,
      Rank   = rank,
      Taxon  = names(counts),
      n_asvs = as.integer(counts)
    )
  })
}) %>%
  mutate(
    Rank   = factor(Rank, levels = tax_ranks, ordered = TRUE),
    Marker = factor(Marker, levels = names(use_ps_list))
  )

# Fraction of ASVs per marker x rank
asv_pct_df <- asv_pct_df %>%
  group_by(Marker, Rank) %>%
  mutate(Pct = n_asvs / sum(n_asvs)) %>%
  ungroup()

# Top N taxa by total ASV count, rest lumped to "Other"
named_taxa <- asv_pct_df %>%
  filter(Taxon != "Unassigned") %>%
  group_by(Taxon) %>%
  summarise(total = sum(n_asvs), .groups = "drop") %>%
  arrange(desc(total))

top_taxa_asv <- head(named_taxa$Taxon, top_n_taxa)

asv_pct_df <- asv_pct_df %>%
  mutate(Taxon_plot = case_when(
    Taxon == "Unassigned" ~ "Unassigned",
    Taxon %in% top_taxa_asv ~ Taxon,
    TRUE ~ "Other"
  ))

asv_pct_df <- asv_pct_df %>%
  group_by(Marker, Rank, Taxon_plot) %>%
  summarise(Pct = sum(Pct), .groups = "drop")

# ASV-count assignment rate table
assign_rate_asv <- asv_pct_df %>%
  filter(Taxon_plot != "Unassigned") %>%
  group_by(Marker, Rank) %>%
  summarise(pct_assigned = round(sum(Pct) * 100, 1), .groups = "drop") %>%
  arrange(Marker, Rank)

save_stats(assign_rate_asv,
           file.path(outdir, "stats", "assignment_rates_asvcount"),
           caption = "% ASVs assigned at each taxonomic rank per marker")