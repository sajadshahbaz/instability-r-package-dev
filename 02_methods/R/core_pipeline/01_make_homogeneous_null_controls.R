#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) cat(sprintf("[%s] ", ts()), sprintf(...), "\n", sep = "")
stop2 <- function(...) {
  cat(sprintf("[FATAL %s] ", ts()), sprintf(...), "\n", sep = "", file = stderr())
  quit(status = 2)
}

read_config_simple <- function(path) {
  if (!file.exists(path)) stop2("Config not found: %s", path)
  x <- readLines(path, warn = FALSE)
  x <- gsub("#.*$", "", x)
  x <- x[nzchar(trimws(x))]
  kv <- list()
  for (ln in x) {
    if (!grepl(":", ln, fixed = TRUE)) next
    if (grepl("^\\s*[A-Za-z0-9_]+\\s*:\\s*$", ln)) next
    m <- regexec("^\\s*([A-Za-z0-9_]+)\\s*:\\s*(.*)\\s*$", ln)
    r <- regmatches(ln, m)[[1]]
    if (length(r) != 3) next
    key <- r[2]
    val <- r[3]
    val <- gsub('^"(.*)"$', "\\1", val)
    kv[[key]] <- val
  }
  kv
}

clean_numeric <- function(x) as.numeric(gsub("[^0-9eE+.-]", "", as.character(x)))

effect_delta <- function(mat, treat_ids, ctrl_ids) {
  mt <- if (length(treat_ids) > 0) rowMeans(mat[, treat_ids, drop = FALSE]) else rep(NA_real_, nrow(mat))
  mc <- if (length(ctrl_ids)  > 0) rowMeans(mat[, ctrl_ids,  drop = FALSE]) else rep(NA_real_, nrow(mat))
  mt - mc
}

sgn <- function(x, eps = 1e-6) {
  if (is.na(x)) return(NA_integer_)
  if (abs(x) <= eps) return(0L)
  if (x > 0) return(1L)
  -1L
}

split_half <- function(ids) {
  ids <- sample(ids)
  n <- length(ids)
  if (n < 2) return(list(a = ids, b = character(0)))
  k <- floor(n / 2)
  list(a = ids[1:k], b = ids[(k + 1):n])
}

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop2("NM2_ROOT is not set")

cfg_path <- file.path(nm2_root, "00_config", "config.yaml")
cfg <- read_config_simple(cfg_path)

vst_path  <- cfg[["vst"]]
meta_path <- cfg[["meta"]]

if (!file.exists(vst_path))  stop2("Missing VST: %s", vst_path)
if (!file.exists(meta_path)) stop2("Missing meta: %s", meta_path)

sid_col  <- cfg[["sample_id_col"]]
cond_col <- cfg[["condition_col"]]
grp_col  <- cfg[["group_col"]]

n_boot <- 200L
seed <- 1L
set.seed(seed)

tables_dir <- file.path(nm2_root, "04_results", "tables")
logs_dir   <- file.path(nm2_root, "04_results", "logs")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir,   showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, "01_homogeneous_null_controls_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

log_msg("Homogeneous null-control generation start")
log_msg("NM2_ROOT = %s", normalizePath(nm2_root))
log_msg("Config   = %s", normalizePath(cfg_path))
log_msg("n_boot=%d seed=%d", n_boot, seed)

meta <- fread(meta_path)
meta[, (sid_col) := trimws(as.character(get(sid_col)))]
meta[, (cond_col) := tolower(trimws(as.character(get(cond_col))))]
meta[, (grp_col)  := tolower(trimws(as.character(get(grp_col))))]

vst_dt <- fread(vst_path)
features <- vst_dt[[1]]
vst_samples <- names(vst_dt)[-1]

common <- intersect(meta[[sid_col]], vst_samples)
if (length(common) == 0) stop2("No sample overlap between meta and VST.")

meta <- meta[get(sid_col) %in% common]
X <- as.matrix(vst_dt[, ..common])
mode(X) <- "numeric"
rownames(X) <- features

run_null_control <- function(cond_label, meta_condition_name, outfile_stub) {
  m <- meta[get(cond_col) == meta_condition_name]
  if (nrow(m) == 0) {
    log_msg("Condition '%s': n=0 (skip)", meta_condition_name)
    return(invisible(NULL))
  }

  treat_all <- m[get(grp_col) == "treatment", get(sid_col)]
  ctrl_all  <- m[get(grp_col) == "control",   get(sid_col)]

  if (length(treat_all) < 4 || length(ctrl_all) < 4) {
    log_msg("Condition '%s': too few treatment/control samples for stable pseudo-splits (treat=%d ctrl=%d) -> skip",
            meta_condition_name, length(treat_all), length(ctrl_all))
    return(invisible(NULL))
  }

  log_msg("COND=%s | pooled: treat=%d ctrl=%d", cond_label, length(treat_all), length(ctrl_all))

  d_pool <- effect_delta(X, treat_all, ctrl_all)
  cond_samples <- m[[sid_col]]
  var_all <- apply(X[, cond_samples, drop = FALSE], 1, var)

  n_feat <- nrow(X)
  hetero_gap_mat <- matrix(NA_real_, nrow = n_feat, ncol = n_boot)
  discord_mat    <- matrix(NA_real_, nrow = n_feat, ncol = n_boot)
  instab_mat     <- matrix(NA_real_, nrow = n_feat, ncol = n_boot)
  signflip_mat   <- matrix(FALSE, nrow = n_feat, ncol = n_boot)

  for (b in seq_len(n_boot)) {
    sp_t <- split_half(treat_all)
    sp_c <- split_half(ctrl_all)

    treat_1 <- sp_t$a
    treat_2 <- sp_t$b
    ctrl_1  <- sp_c$a
    ctrl_2  <- sp_c$b

    # require both pseudo-groups to have both treatment and control
    if (length(treat_1) == 0 || length(treat_2) == 0 || length(ctrl_1) == 0 || length(ctrl_2) == 0) next

    d1 <- effect_delta(X, treat_1, ctrl_1)
    d2 <- effect_delta(X, treat_2, ctrl_2)

    n1 <- length(treat_1) + length(ctrl_1)
    n2 <- length(treat_2) + length(ctrl_2)
    w1 <- n1 / (n1 + n2)
    w2 <- n2 / (n1 + n2)
    d_strat <- w1 * d1 + w2 * d2

    sgn1 <- vapply(d1, sgn, integer(1))
    sgn2 <- vapply(d2, sgn, integer(1))

    sign_flip   <- (sgn1 != 0L) & (sgn2 != 0L) & (sgn1 != sgn2)
    hetero_gap  <- abs(d1 - d2)
    discordance <- abs(d_pool - d_strat)
    instability <- hetero_gap + discordance

    hetero_gap_mat[, b] <- hetero_gap
    discord_mat[, b]    <- discordance
    instab_mat[, b]     <- instability
    signflip_mat[, b]   <- sign_flip
  }

  out <- data.table(
    gene = rownames(X),
    condition = cond_label,
    delta_pool = d_pool,
    hetero_gap = rowMeans(hetero_gap_mat, na.rm = TRUE),
    discordance = rowMeans(discord_mat, na.rm = TRUE),
    instability_score = rowMeans(instab_mat, na.rm = TRUE),
    sign_flip_rate = rowMeans(signflip_mat, na.rm = TRUE),
    var_proxy = var_all,
    n_boot = n_boot
  )

  out_path <- file.path(tables_dir, paste0(outfile_stub, "_null_control_all.tsv"))
  sum_path <- file.path(tables_dir, paste0(outfile_stub, "_null_control_summary.tsv"))

  fwrite(out, out_path, sep = "\t")

  summary <- out[, .(
    n_features = .N,
    median_hetero_gap = median(hetero_gap, na.rm = TRUE),
    median_discordance = median(discordance, na.rm = TRUE),
    median_instability = median(instability_score, na.rm = TRUE),
    p90_instability = quantile(instability_score, 0.90, na.rm = TRUE),
    p95_instability = quantile(instability_score, 0.95, na.rm = TRUE),
    mean_sign_flip_rate = mean(sign_flip_rate, na.rm = TRUE),
    cor_instability_vs_var = suppressWarnings(cor(instability_score, var_proxy, use = "pairwise.complete.obs"))
  )]
  fwrite(summary, sum_path, sep = "\t")

  log_msg("Saved: %s", out_path)
  log_msg("Saved: %s", sum_path)
}

run_null_control("HT",  "ht",      "ht")
run_null_control("LT",  "lt",      "lt")
run_null_control("OSM", "osmotic", "osmotic")

log_msg("Homogeneous null-control generation done.")
