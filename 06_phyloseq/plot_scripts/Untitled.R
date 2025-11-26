# Required packages
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
# optional for a quick pairs plot
if (!requireNamespace("GGally", quietly = TRUE)) install.packages("GGally")
library(GGally)

# 1) helper: count unique taxa per sample at a given rank
count_taxa_per_sample <- function(ps, rank = "Genus", min_abundance = 1,
                                  include_unclassified = FALSE) {
  # ps: phyloseq object
  # rank: string, e.g. "Genus" or "Species"
  # min_abundance: minimum reads to consider taxon present in sample
  # include_unclassified: TRUE treat NA/"" as "Unclassified", FALSE drop them
  df <- psmelt(ps) %>%
    filter(Abundance >= min_abundance) %>%
    mutate(taxon = as.character(.data[[rank]])) %>%
    mutate(taxon = ifelse(is.na(taxon) | taxon == "",
                          ifelse(include_unclassified, "Unclassified", NA_character_),
                          taxon)) %>%
    filter(!is.na(taxon)) %>%
    group_by(Sample) %>%
    summarise(n_taxa = n_distinct(taxon), .groups = "drop")
  
  # ensure every sample in ps is present (0 if none)
  all_samples <- sample_names(ps)
  out <- tibble::tibble(Sample = all_samples) %>%
    left_join(df, by = "Sample") %>%
    mutate(n_taxa = ifelse(is.na(n_taxa), 0L, as.integer(n_taxa)))
  
  out
}

# 2) assemble counts from a named list of phyloseq objects
compute_counts_matrix <- function(ps_list, rank = "Genus", min_abundance = 1,
                                  include_unclassified = FALSE) {
  # ps_list: named list, e.g. list(ps1 = ps1, ps2 = ps2, ps3 = ps3, ps4 = ps4)
  if (is.null(names(ps_list)) || any(names(ps_list) == "")) {
    stop("ps_list must be a named list (names will become column names).")
  }
  
  # take sample set from first object; ensure others contain same sample names
  samples_ref <- sample_names(ps_list[[1]])
  # optional: check that all sample sets match (warning if not)
  for (nm in names(ps_list)[-1]) {
    if (!setequal(samples_ref, sample_names(ps_list[[nm]]))) {
      warning("Sample names differ between first object and ", nm,
              ". Intersection will be used.")
    }
  }
  
  # compute counts per object
  counts_df <- tibble::tibble(Sample = samples_ref)
  for (nm in names(ps_list)) {
    temp <- count_taxa_per_sample(ps_list[[nm]], rank, min_abundance, include_unclassified)
    # align rows by Sample (use intersection)
    temp <- temp %>% filter(Sample %in% samples_ref)
    counts_df <- counts_df %>% left_join(temp, by = "Sample")
    colnames(counts_df)[ncol(counts_df)] <- nm
  }
  # ensure integer columns
  for (nm in names(ps_list)) counts_df[[nm]] <- as.integer(counts_df[[nm]])
  counts_df
}

# 3) Quick pairs plot (GGally)
plot_counts_pairs_ggpairs <- function(counts_df, log_transform = FALSE) {
  # counts_df: output of compute_counts_matrix (Sample column + object columns)
  mat <- counts_df %>% select(-Sample)
  if (log_transform) mat <- log10(mat + 1)
  # GGally will show histograms on diagonal, scatter+smooth and correlations
  p <- GGally::ggpairs(mat,
                       upper = list(continuous = wrap("cor", size = 4)),
                       lower = list(continuous = wrap("smooth", alpha = 0.3, se = FALSE)),
                       diag = list(continuous = "barDiag"))
  p
}

# 4) Custom pairwise ggplot with lm + R^2 annotation (saves list of plots)
pairwise_scatter_plots <- function(counts_df, log_transform = FALSE, prefix = "pairplot") {
  # returns a named list of ggplot objects
  mat <- counts_df %>% select(-Sample)
  if (log_transform) {
    mat <- mat %>% mutate_all(~log10(.x + 1))
  }
  
  obj_names <- colnames(mat)
  combos <- combn(obj_names, 2, simplify = FALSE)
  
  plots <- list()
  for (cmb in combos) {
    a <- cmb[1]; b <- cmb[2]
    dat <- tibble::tibble(x = mat[[a]], y = mat[[b]], Sample = counts_df$Sample)
    fit <- lm(y ~ x, data = dat)
    r2 <- summary(fit)$r.squared
    p <- ggplot(dat, aes(x = x, y = y)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
      labs(x = a, y = b, title = paste(a, "vs", b),
           subtitle = paste0("R² = ", formatC(r2, digits = 3))) +
      theme_minimal()
    name <- paste(a, "_vs_", b, sep = "")
    plots[[name]] <- p
  }
  plots
}

########################
# Example usage:


Leray_ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/Leray_dadaOriginal.RDS")
Folmer_ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/Folmer_dadaOriginal.RDS")
ps_18S <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/18S_dadaOriginal.RDS")
morph_ps <- readRDS("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/ps_objects/morph_ps_integer.RDS")


morph_names <- sample_names(morph_ps)
comp_folmer_ps <- subset_samples(Folmer_ps, Station_ID %in% morph_names)
comp_leray_ps <- subset_samples(Leray_ps, Station_ID %in% morph_names)
comp_18S_ps <- subset_samples(ps_18S, Station_ID %in% morph_names)

comp_folmer_ps <- merge_samples(comp_folmer_ps, "Station_ID")
comp_leray_ps <- merge_samples(comp_leray_ps, "Station_ID")
comp_18S_ps <- merge_samples(comp_18S_ps, "Station_ID")


ps_list <- list(folmer = comp_folmer_ps,
                leray = comp_leray_ps,
                SSU = comp_18S_ps,
                morph = morph_ps,
                leray_folm = uhhhh)


# Taxa count function

taxa_sum_per_sample <- function(ps, tax_rank = "Genus"){
  hello <- ps %>% 
    tax_glom(taxrank = tax_rank) %>% 
    transform_sample_counts(fun = function(x) (x > 0) * 1L) %>% 
    otu_table() %>% 
    as.matrix() -> mat
  
  
  if (taxa_are_rows(otu_table(ps))) {
    # number of samples in which each Class is present
    hello <- colSums(mat)
  } else {
    hello <- rowSums(mat)
  }
  hello <- as.data.frame(hello)
  colnames(hello) <- c(tax_rank)
  return(hello)
}

leray_class <- (taxa_sum_per_sample(comp_leray_ps, "Class"))

# Loop to count all taxa levels for each marker
outputs <- vector("list", length(ps_list))   # pre-allocate
names(outputs) <- names(ps_list)
for(i in seq_along(ps_list)){
  nm <- names(ps_list)[i]
  ps <- ps_list[[i]]
  
  taxa_count <- data.frame(row.names = sample_names(ps))
  ranks <- rank_names(ps)
  for(rank in ranks){
    countcol <- taxa_sum_per_sample(ps, rank)
    df_to_bind <- countcol[rownames(taxa_count), , drop = FALSE] 
    colnames(df_to_bind) <- rank
    taxa_count <- cbind(taxa_count, df_to_bind)
    outputs[[nm]] <- taxa_count
    #  taxa_count <- merge(taxa_count, countcol, by = "row.names", all.x = TRUE)
  }
}

str(outputs)




a <- outputs[[1]]
b <- outputs[[4]]

# basic checks & extraction (robust to factor/character)
if(!"Class" %in% colnames(a) || !"Class" %in% colnames(b)) {
  stop("No 'Class' column found in one of the data.frames.")
}

# coerce to numeric & align by sample names
vec_a <- as.numeric(as.character(a[rownames(a), "Class"]))
vec_b <- as.numeric(as.character(b[rownames(b), "Class"]))

common <- intersect(rownames(a), rownames(b))

plot_df <- data.frame(
  Sample = common,
  a_Class = as.numeric(as.character(a[common, "Species"])),
  b_Class = as.numeric(as.character(b[common, "Species"])),
  stringsAsFactors = FALSE,
  row.names = NULL
)


# remove NAs (if any)
plot_df <- na.omit(plot_df)

# fit linear model and get stats
fit <- lm(b_Class ~ a_Class, data = plot_df)
s <- summary(fit)
r2 <- s$r.squared
pval <- coef(s)[2,4]    # p-value for slope

# basic scatter + lm line + annotation
p <- ggplot(plot_df, aes(x = a_Class, y = b_Class)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  labs(
    x = "Unique Classes (outputs[[1]])",
    y = "Unique Classes (outputs[[2]])",
    title = "Per-sample unique Classes: outputs[[1]] vs outputs[[2]]"
  ) +
  theme_minimal() +
  #annotate R^2 and p-value at top-right of panel
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("R² = ", formatC(r2, digits = 3), "\n", "p = ", format.pval(pval, digits = 3)),
           hjust = 1.1, vjust = 1.1, size = 3.5)

print(p)









library(ggplot2)

plot_compare_safe_fixed <- function(df1, df2, rank = "Class",
                                    df1_name = NULL, df2_name = NULL,
                                    label_samples = FALSE, log10 = FALSE,
                                    min_points_for_lm = 2) {
  # df1, df2: data.frames with rownames = sample IDs and a column named `rank`
  # rank: column name to compare (character)
  # df1_name / df2_name: optional labels for axes/title (defaults to "df1"/"df2")
  # label_samples: use ggrepel to label points if TRUE and ggrepel is installed
  # log10: plot log10(x+1), log10(y+1) instead of raw counts
  # min_points_for_lm: minimum number of non-NA pairs required to fit lm
  
  if (is.null(df1_name)) df1_name <- deparse(substitute(df1))
  if (is.null(df2_name)) df2_name <- deparse(substitute(df2))
  
  # sanity checks
  if (!rank %in% colnames(df1)) stop("rank '", rank, "' not found in df1")
  if (!rank %in% colnames(df2)) stop("rank '", rank, "' not found in df2")
  if (is.null(rownames(df1)) || is.null(rownames(df2))) {
    stop("Both data.frames must have rownames as sample IDs")
  }
  
  # get intersection of samples and preserve order in df1 (but use intersect to ensure common)
  common <- intersect(rownames(df1), rownames(df2))
  if (length(common) == 0) {
    return(ggplot() + ggtitle("No common samples between the two data.frames") + theme_minimal())
  }
  
  # extract safely by row index (this avoids name-lookup issues)
  x_raw <- df1[common, rank, drop = TRUE]
  y_raw <- df2[common, rank, drop = TRUE]
  
  # safe numeric coercion function
  safe_num <- function(v) {
    if (is.factor(v)) v <- as.character(v)
    # convert empty strings to NA first
    v[v == ""] <- NA
    v_num <- suppressWarnings(as.numeric(as.character(v)))
    return(v_num)
  }
  
  x <- safe_num(x_raw)
  y <- safe_num(y_raw)
  
  # make data.frame and drop pairs where either is NA
  df <- data.frame(Sample = common, x = x, y = y, stringsAsFactors = FALSE)
  df <- df[!is.na(df$x) & !is.na(df$y), , drop = FALSE]
  
  if (nrow(df) == 0) {
    return(ggplot() + ggtitle("No non-NA paired values to plot") + theme_minimal())
  }
  
  # optional log transform
  if (log10) {
    df$x <- log10(df$x + 1)
    df$y <- log10(df$y + 1)
  }
  
  # base scatter
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(size = 2) +
    labs(
      x = paste0("Unique ", rank, " (", df1_name, ")"),
      y = paste0("Unique ", rank, " (", df2_name, ")"),
      title = paste0("Per-sample unique ", rank, ": ", df1_name, " vs ", df2_name)
    ) +
    theme_minimal()
  
  # fit lm only if we have enough points
  if (nrow(df) >= min_points_for_lm) {
    fit <- lm(y ~ x, data = df)
    s <- summary(fit)
    r2 <- s$r.squared
    pval <- coef(s)[2, 4]
    # add regression line and annotation
    p <- p +
      geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
      annotate("text", x = Inf, y = Inf,
               label = paste0("n=", nrow(df), "  R²=", formatC(r2, digits = 3),
                              "\n", "p=", format.pval(pval, digits = 3),
                              "\n", "slope=", formatC(coef(fit)[2], digits = 3)),
               hjust = 1.05, vjust = 1.05, size = 3.5)
  } else {
    p <- p + ggtitle(paste0("Only ", nrow(df), " non-NA point(s); no linear fit"))
  }
  
  # optional point labels using ggrepel (if requested and available)
  if (label_samples) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = df, aes(label = Sample), size = 3)
    } else {
      warning("ggrepel not installed; install.packages('ggrepel') to use non-overlapping labels")
      p <- p + geom_text(aes(label = Sample), vjust = -0.5, size = 2.5)
    }
  }
  
  p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey40")
  p <- p + theme_bw()
  
  return(p)
}

# outputs is a list of data.frames; each df has rownames = sample names and column "Class"
p <- plot_compare_safe_fixed(outputs[[2]], outputs[[1]],
                             rank = "Species",
                             df1_name = "Leray", df2_name = "Folmer",
                             label_samples = TRUE, log10 = FALSE)
#p1 <- p + annotate("text", x = Inf, y = Inf, 
 #                 label = paste0("slope = ", formatC(coef(fit)[2], digits = 3)),
  #                hjust = 0.001, vjust = 0.001, size = 3.5)

print(p)

outdir = "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/06_phyloseq/plots/"
ggsave(paste0(outdir,"Morph_vs_Leray+Folmer_Species.png"), plot = p,
       width = 1000, height = 800, device = png, units = "px", scale = 3)







str(outputs)


a <- get_taxa_unique(comp_folmer_ps, taxonomic.rank = "Species")
b <- get_taxa_unique(comp_leray_ps, taxonomic.rank = "Species")
c <- c(a,b)
length(a)
length(b)
length(unique(c))


comp_leray_glom_ps <- tax_glom(comp_leray_ps, taxrank = "Species")
comp_folmer_glom_ps <- tax_glom(comp_folmer_ps, taxrank = "Species")
uhhhh <- merge_phyloseq(comp_leray_glom_ps, comp_folmer_glom_ps)
c1 <- get_taxa_unique(uhhhh, taxonomic.rank = "Species" )
length(c1)


hello <- comp_leray_ps %>% 
  tax_glom(taxrank = "Class") %>% 
  transform_sample_counts(fun = function(x) (x > 0) * 1L) %>% 
  otu_table() %>% 
  as.matrix() -> mat


if (taxa_are_rows(otu_table(comp_leray_ps))) {
  # number of samples in which each Class is present
  hello <- colSums(mat)
} else {
  hello <- rowSums(mat)
}
hello 


(otu_table(transform_sample_counts(tax_glom(comp_leray_ps, taxrank = "Class"))))




count_taxa_per_sample(comp_leray_ps, rank = "Class")






ps_list <- list(folmer = comp_folmer_ps, leray = comp_leray_ps, 
                SSU = comp_18S_ps, morph = morph_ps)
counts_df <- compute_counts_matrix(ps_list, rank = "Genus", min_abundance = 1, include_unclassified = FALSE)
View(counts_df)  # Sample | obj1 | obj2 | obj3 | obj4

# Quick matrix plot:
p_pairs <- plot_counts_pairs_ggpairs(counts_df, log_transform = FALSE)
print(p_pairs)

# Or get individual pair plots:
pair_plots <- pairwise_scatter_plots(counts_df, log_transform = TRUE)
# show one:
print(pair_plots[[1]])
# save all:
for (nm in names(pair_plots)) ggsave(filename = paste0("plots/", nm, ".png"), plot = pair_plots[[nm]], width = 5, height = 5)
########################
