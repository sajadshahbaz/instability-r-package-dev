#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

# ---------------------------
# utilities
# ---------------------------

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
  cat(sprintf("[%s] ", ts()), sprintf(...), "\n", sep = "")
}

stop2 <- function(...) {
  cat(sprintf("[FATAL %s] ", ts()), sprintf(...), "\n", sep = "", file = stderr())
  quit(status = 2)
}

assert_readable <- function(path, label) {
  if (!file.exists(path))
    stop2("%s path does not exist: %s", label, path)
  if (file.access(path, 4) != 0)
    stop2("%s path not readable: %s", label, path)
}



# ---------------------------
# config
# ---------------------------

config_path <- "00_config/config.yaml"

if (!file.exists(config_path))
  stop2("Missing config.yaml: %s", config_path)

cfg <- tryCatch(
  yaml::read_yaml(config_path),
  error = function(e) stop2("Failed to parse YAML: %s", e$message)
)

if (is.null(cfg$paths))
  stop2("Config missing 'paths' section")

paths <- cfg$paths

vst_path    <- paths$vst
meta_path   <- paths$meta
counts_path <- paths$counts

if (is.null(vst_path) || is.null(meta_path))
  stop2("Config must define paths$vst and paths$meta")

log_msg("NM2 IO start")
log_msg("Config: %s", normalizePath(config_path))
log_msg("VST:   %s", vst_path)
log_msg("META:  %s", meta_path)

if (!is.null(counts_path) && nzchar(counts_path)) {
  log_msg("COUNTS(optional): %s", counts_path)
}

# ---------------------------
# validate input paths
# ---------------------------

assert_readable(vst_path,  "VST")
assert_readable(meta_path, "META")

if (!is.null(counts_path) && nzchar(counts_path)) {
  if (!file.exists(counts_path) || file.access(counts_path, 4) != 0) {
    stop2("Counts path provided but not readable: %s", counts_path)
  }
}

# ---------------------------
# read files
# ---------------------------

meta <- tryCatch(
  fread(meta_path),
  error = function(e) stop2("Failed to read META: %s", e$message)
)

vst <- tryCatch(
  fread(vst_path),
  error = function(e) stop2("Failed to read VST: %s", e$message)
)

if (ncol(vst) < 2)
  stop2("VST file malformed (ncol < 2): %s", vst_path)

if (!"sample_id" %in% names(meta))
  stop2("META missing required column: sample_id")

# ---------------------------
# sample alignment
# ---------------------------

vst_samples  <- colnames(vst)[-1]
meta_samples <- meta$sample_id

overlap <- intersect(vst_samples, meta_samples)

log_msg("Meta samples: %d", length(meta_samples))
log_msg("VST samples:  %d", length(vst_samples))
log_msg("Overlap:      %d", length(overlap))

meta_only <- setdiff(meta_samples, vst_samples)
vst_only  <- setdiff(vst_samples, meta_samples)

log_msg("Meta-only samples: %d", length(meta_only))
log_msg("VST-only samples:  %d", length(vst_only))

if (length(overlap) == 0)
  stop2("No overlapping sample IDs between metadata and VST columns")

# ---------------------------
# per-condition summary
# ---------------------------

if (!"condition" %in% names(meta))
  stop2("META missing required column: condition")

if (!"baseline_block" %in% names(meta))
  stop2("META missing required column: baseline_block")

if (!"group_type" %in% names(meta))
  stop2("META missing required column: group_type")

meta[, condition      := tolower(trimws(condition))]
meta[, baseline_block := tolower(trimws(baseline_block))]
meta[, group_type     := tolower(trimws(group_type))]

conds <- unique(meta$condition)

summary_lines <- list()

for (cond in conds) {
  msub <- meta[condition == cond]
  
  log_msg("---- Condition: %s  (n=%d)", cond, nrow(msub))
  
  by_block <- msub[, .N, by = baseline_block]
  by_block_group <- msub[, .N, by = .(baseline_block, group_type)]
  
  print(by_block)
  print(by_block_group)
  
  summary_lines[[cond]] <- list(
    n_total = nrow(msub),
    by_block = by_block,
    by_block_group = by_block_group
  )
}

# ---------------------------
# write summary log
# ---------------------------

out_dir <- "04_results/logs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(out_dir, "00_io_summary.txt")
sink(summary_file)
cat("NM2 IO Summary\n")
cat("Generated:", ts(), "\n\n")

for (cond in names(summary_lines)) {
  cat("Condition:", cond, "\n")
  cat("Total samples:", summary_lines[[cond]]$n_total, "\n")
  cat("\nBy baseline_block:\n")
  print(summary_lines[[cond]]$by_block)
  cat("\nBy baseline_block x group_type:\n")
  print(summary_lines[[cond]]$by_block_group)
  cat("\n---------------------------\n\n")
}

sink()

log_msg("NM2 IO done. Summary written to: %s", summary_file)
