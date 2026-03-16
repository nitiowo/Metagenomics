# run_all.R

# edit config.R to change paths and output_root

source("config.R")

source("setup.R")

script_steps <- c(
  exploratory = "exploratory.R",
  composition = "composition.R",
  alpha = "alpha.R",
  beta = "beta.R",
  differential = "differential.R",
  heatmaps = "heatmaps.R",
  overlap = "overlap.R",
  geographic = "geographic.R",
  focal_taxon = "focal_taxon.R",
  varpart = "varpart.R",
  trebitz_compare = "trebitz_compare.R"


# Set step to FALSE to skip
run_steps <- setNames(rep(TRUE, length(script_steps)), names(script_steps))


for (nm in names(run_steps)) {
  cat("  ", nm, "->", ifelse(isTRUE(run_steps[[nm]]), "RUN", "SKIP"), "\n")
}

for (nm in names(script_steps)) {
  s <- script_steps[[nm]]

  if (!isTRUE(run_steps[[nm]])) {
    cat("\n--- skipping:", nm, "---\n")
    next
  }

  if (file.exists(s)) {
    cat("\n--- running:", nm, "(", s, ") ---\n")
    tryCatch(source(s), error = function(e) cat("  ERROR:", e$message, "\n"))
  } else {
    cat("  skipping (not found):", s, "\n")
  }
}

cat("  done. output in:", output_root, "\n")
