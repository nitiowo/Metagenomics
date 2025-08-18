

manifest <- read_tsv(manifest_path, col_types = cols(
  sample = col_character(),
  reverse_bin = col_character(),
  out_R1_path = col_character(),
  out_R2_path = col_character(),
  CA_exit_status = col_logical()
))

# Helper: safely run a shell command and return trimmed stdout (or NA on error)
run_shell <- function(cmd) {
  out <- tryCatch({
    res <- system2("bash", c("-lc", cmd), stdout = TRUE, stderr = FALSE)
    if (length(res) == 0) return(NA_character_)
    paste(res, collapse = "\n")
  }, error = function(e) NA_character_)
  if (is.na(out)) return(NA_character_)
  trimws(out)
}

# # Count reads by calling: if file exists -> zcat 'file' | wc -l -> /4
# count_reads_fastq <- function(fq) {
#   if (is.na(fq) || fq == "" || !file.exists(fq)) return(0L)
#   cmd <- sprintf("zcat %s | wc -l", shQuote(fq))
#   lines_s <- run_shell(cmd)
#   if (is.na(lines_s) || lines_s == "") return(0L)
#   lines_n <- as.numeric(lines_s)
#   if (is.na(lines_n)) return(0L)
#   reads <- as.integer(lines_n / 4)
#   return(reads)
# }

# Average read length using awk on sequence lines (every 4th line, NR%4==2)
avg_len_fastq <- function(fq) {
  if (is.na(fq) || fq == "" || !file.exists(fq)) return(0)
  # awk returns a single number (or 0)
  awk_cmd <- "awk 'NR % 4 == 2 { total += length($0); count++ } END { if (count>0) printf(\"%.6f\", total/count); else print 0 }'"
  cmd <- sprintf("zcat %s | %s", shQuote(fq), awk_cmd)
  res <- run_shell(cmd)
  if (is.na(res) || res == "") return(0)
  num <- as.numeric(res)
  if (is.na(num)) return(0)
  return(num)
}

# -------------------------
# Compute counts & avg lengths for each row
# -------------------------
message("Computing counts and average lengths for manifest rows (this may take a few minutes)...")

# To avoid re-running zcat many times in large manifests, you could parallelize this loop.
rows <- nrow(manifest)
results <- manifest %>%
  select(-(out_R1_path:CA_exit_status)) %>% 
  mutate(
    Count_R = 0L,
    Count_F = 0L,
    AvgLen_R = 0.0,
    AvgLen_F = 0.0
  )

for (i in seq_len(rows)) {
  samp <- manifest$sample[1]
  primer <- manifest$reverse_bin[1]
  r1 <- manifest$out_R1_path[1]   # forward read path in your manifest (R1)
  r2 <- manifest$out_R2_path[1]   # reverse read path (R2)
  # Note: by convention in this script Count_R corresponds to R2_path, Count_F to R1_path
  cr <- as.integer(countLines(r2))
  cf <- as.integer(countLines(r2))
  ar <- avg_len_fastq(r2)
  af <- avg_len_fastq(r1)
  results$Count_R[i] <- cr
  results$Count_F[i] <- cf
  results$AvgLen_R[i] <- ar
  results$AvgLen_F[i] <- af
  if (i %% 10 == 0 || i == rows) {
    message(sprintf("  processed %d / %d rows", i, rows))
  }
}

# -------------------------
# Per-sample FR_countsMatch: for each sample, is Count_R == Count_F for every primer?
# -------------------------
fr_flag <- results %>%
  group_by(Sample) %>%
  summarise(FR_CountsMatch = ifelse(all(Count_R == Count_F), "YES", "NO")) 

# -------------------------
# Tot_countsMatch: compare totals to Raw_R paths if provided
# If manifest includes Raw_R1_path/Raw_R2_path for each sample (or at least one per sample), compute raw counts and compare.
# -------------------------
# find one Raw path per sample (if present)
raw_paths <- manifest %>%
  group_by(Sample) %>%
  summarise(
    Raw_R1_path = first(na.omit(Raw_R1_path)),
    Raw_R2_path = first(na.omit(Raw_R2_path))
  )

# compute totals per sample from demuxed results
totals <- results %>%
  group_by(Sample) %>%
  summarise(Tot_Count_R = sum(Count_R), Tot_Count_F = sum(Count_F))

# join raw paths and compute raw counts if paths exist
raw_counts <- raw_paths %>%
  left_join(totals, by = "Sample") %>%
  rowwise() %>%
  mutate(
    Raw_Count_R = if (!is.na(Raw_R2_path) && Raw_R2_path != "" && file.exists(Raw_R2_path)) count_reads_fastq(Raw_R2_path) else NA_integer_,
    Raw_Count_F = if (!is.na(Raw_R1_path) && Raw_R1_path != "" && file.exists(Raw_R1_path)) count_reads_fastq(Raw_R1_path) else NA_integer_,
    Tot_countsMatch = if (!is.na(Raw_Count_R) && !is.na(Raw_Count_F)) {
      if (Raw_Count_R == Tot_Count_R && Raw_Count_F == Tot_Count_F) "YES" else "NO"
    } else {
      NA_character_
    }
  ) %>%
  ungroup() %>%
  select(Sample, Tot_countsMatch)

# -------------------------
# Combine flags and prepare wide summaries
# -------------------------
summary_full <- results %>%
  mutate(FR_countsMatch = ifelse(Count_R == Count_F, "YES", "NO")) %>%
  left_join(fr_flag, by = "Sample", suffix = c("", ".sample")) %>%
  left_join(raw_counts, by = "Sample")

# create a short name mapping for the known primers (customize if you have other primer names)
short_map <- c("HCO2198" = "HCO", "SSUR22" = "SSU", "unknown" = "unk")
summary_full <- summary_full %>%
  mutate(Short = ifelse(Primer %in% names(short_map), short_map[Primer], Primer))

# counts summary: pivot to have Count columns per primer (Count_HCO etc.) using Count_R
counts_wide <- summary_full %>%
  select(Sample, Short, Count_R) %>%
  mutate(Short = paste0("Count_", Short)) %>%
  pivot_wider(names_from = Short, values_from = Count_R, values_fill = 0) %>%
  left_join(
    summary_full %>%
      group_by(Sample) %>%
      summarise(FR_CountsMatch = first(FR_CountsMatch), Tot_countsMatch = first(Tot_countsMatch)),
    by = "Sample"
  )

# read length summary: pivot AvgLen_R and AvgLen_F into columns like readLen_HCO_R, readLen_HCO_F
len_wide_R <- summary_full %>%
  select(Sample, Short, AvgLen_R) %>%
  mutate(Short = paste0("readLen_", Short, "_R")) %>%
  pivot_wider(names_from = Short, values_from = AvgLen_R, values_fill = 0)

len_wide_F <- summary_full %>%
  select(Sample, Short, AvgLen_F) %>%
  mutate(Short = paste0("readLen_", Short, "_F")) %>%
  pivot_wider(names_from = Short, values_from = AvgLen_F, values_fill = 0)

readlen_wide <- len_wide_R %>%
  full_join(len_wide_F, by = "Sample") %>%
  arrange(Sample)

# Re-order columns to the conventional order if HCO, SSU, unk exist
desired_count_cols <- c("Sample",
                        paste0("Count_", c("HCO","SSU","unk")),
                        "FR_CountsMatch","Tot_countsMatch")

# make sure columns exist (if a primer missing, fill zeros)
for (col in desired_count_cols) {
  if (!col %in% colnames(counts_wide)) counts_wide[[col]] <- 0
}
counts_wide <- counts_wide %>%
  select(all_of(desired_count_cols))

# For readlen, build desired column order
desired_len_cols <- c("Sample")
for (s in c("HCO","SSU","unk")) {
  desired_len_cols <- c(desired_len_cols, paste0("readLen_", s, "_R"), paste0("readLen_", s, "_F"))
}
for (col in desired_len_cols) {
  if (!col %in% colnames(readlen_wide)) readlen_wide[[col]] <- 0
}
readlen_wide <- readlen_wide %>% select(all_of(desired_len_cols))

# -------------------------
# Write TSV outputs
# -------------------------
count_out <- file.path(outdir, "rev_count_summary.tsv")
len_out <- file.path(outdir, "rev_readLen_summary.tsv")

write_tsv(counts_wide, count_out)
write_tsv(readlen_wide, len_out)

message("Wrote: ", count_out)
message("Wrote: ", len_out)
