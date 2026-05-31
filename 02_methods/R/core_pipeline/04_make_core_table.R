#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- Sys.getenv("NM2_ROOT")
if (ROOT == "") {
  ROOT <- normalizePath(getwd(), mustWork = TRUE)
  Sys.setenv(NM2_ROOT = ROOT)
}

out_tab <- file.path(ROOT, "04_results", "tables")
out_log <- file.path(ROOT, "04_results", "logs")

dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log, recursive = TRUE, showWarnings = FALSE)

stop2 <- function(...) {
  cat(sprintf("[FATAL] "), sprintf(...), "\n", sep = "", file = stderr())
  quit(status = 2)
}

tsv_ok <- function(p) {
  if (!file.exists(p)) stop2("Missing: %s", p)
}

in_metrics  <- file.path(out_tab, "pool_vs_strat_metrics.tsv")
in_headline <- file.path(out_tab, "tail_headline_summary.tsv")

tsv_ok(in_metrics)
tsv_ok(in_headline)

m <- fread(in_metrics)
h <- fread(in_headline)

x <- merge(m, h, by = "condition", all = TRUE)

keep <- c(
  "condition",
  "n_features",
  "spearman_instability_vs_effect_divergence",
  "spearman_pool_vs_strat",
  "median_abs_pool_minus_strat",
  "q95_abs_pool_minus_strat",
  "sign_flip_rate",
  "j_gap_vs_dis_K200",
  "j_gap_vs_dis_K1000",
  "flip_top200",
  "flip_top1000"
)

miss <- setdiff(keep, names(x))
if (length(miss)) {
  stop2("Core table missing cols: %s", paste(miss, collapse = ", "))
}

x <- x[, ..keep]
setorder(x, condition)

out <- file.path(out_tab, "core_instability_signal_table.tsv")

fwrite(x, out, sep = "\t")

cat("[OK] Saved:", out, "\n")
