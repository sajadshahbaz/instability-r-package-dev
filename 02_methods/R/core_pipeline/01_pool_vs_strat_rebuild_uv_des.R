#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) cat(sprintf("[%s] ", ts()), sprintf(...), "\n", sep="")
stop2 <- function(...) { cat(sprintf("[FATAL %s] ", ts()), sprintf(...), "\n", sep="", file=stderr()); quit(status=2) }

# --- tiny YAML reader (same style as 00_io.R) ---
read_config_simple <- function(path) {
  if (!file.exists(path)) stop2("Config not found: %s", path)
  x <- readLines(path, warn=FALSE)
  x <- gsub("#.*$", "", x)
  x <- x[nzchar(trimws(x))]
  kv <- list()
  for (ln in x) {
    if (!grepl(":", ln, fixed=TRUE)) next
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

# --- helpers ---
safe_mean <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm=TRUE)
sgn <- function(x, eps=1e-6) {
  if (is.na(x)) return(NA_integer_)
  if (abs(x) <= eps) return(0L)
  if (x > 0) return(1L)
  -1L
}

effect_delta <- function(mat, treat_ids, ctrl_ids) {
  # mat: features x samples numeric matrix (or data.table subset)
  mt <- if (length(treat_ids) > 0) rowMeans(mat[, treat_ids, drop=FALSE]) else rep(NA_real_, nrow(mat))
  mc <- if (length(ctrl_ids)  > 0) rowMeans(mat[, ctrl_ids,  drop=FALSE]) else rep(NA_real_, nrow(mat))
  mt - mc
}

# --- main ---
nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop2("NM2_ROOT is not set")
cfg_path <- file.path(nm2_root, "00_config", "config.yaml")
cfg <- read_config_simple(cfg_path)

vst_path <- cfg[["vst"]]
meta_path <- cfg[["meta"]]
if (!file.exists(vst_path)) stop2("Missing VST: %s", vst_path)
if (!file.exists(meta_path)) stop2("Missing meta: %s", meta_path)

sid_col  <- cfg[["sample_id_col"]]
cond_col <- cfg[["condition_col"]]
sub_col  <- cfg[["subgroup_col"]]
grp_col  <- cfg[["group_col"]]

boot <- as.integer(cfg[["boot"]]); if (is.na(boot)) boot <- 500L
seed <- as.integer(cfg[["seed"]]); if (is.na(seed)) seed <- 1L
topn <- as.integer(cfg[["topn"]]); if (is.na(topn)) topn <- 150L

tables_dir <- file.path(nm2_root, "04_results", "tables")
logs_dir   <- file.path(nm2_root, "04_results", "logs")

dir.create(tables_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(logs_dir, showWarnings=FALSE, recursive=TRUE)

log_file <- file.path(logs_dir, "01_pool_vs_strat_log.txt")
sink(log_file, split=TRUE); on.exit(sink(), add=TRUE)

log_msg("NM2 pool-vs-strat start")
log_msg("Config: %s", normalizePath(cfg_path))
log_msg("Seed=%d Boot=%d TopN=%d", seed, boot, topn)

meta <- fread(meta_path)
meta[, (sid_col) := trimws(as.character(get(sid_col)))]
meta[, (cond_col) := tolower(trimws(as.character(get(cond_col)))) ]
meta[, (sub_col)  := tolower(trimws(as.character(get(sub_col)))) ]
meta[, (grp_col)  := tolower(trimws(as.character(get(grp_col)))) ]

vst_dt <- fread(vst_path)
feat_col <- names(vst_dt)[1]
features <- vst_dt[[1]]
vst_samples <- names(vst_dt)[-1]

# align samples
common <- intersect(meta[[sid_col]], vst_samples)
if (length(common) == 0) stop2("No sample overlap between meta and VST.")
meta <- meta[get(sid_col) %in% common]
setkeyv(meta, sid_col)

# build matrix features x samples (numeric)
X <- as.matrix(vst_dt[, ..common])
mode(X) <- "numeric"
rownames(X) <- features

run_condition <- function(cond) {
  m <- meta[get(cond_col) == cond]
  if (nrow(m) == 0) {
    log_msg("Condition '%s': n=0 (skip)", cond)
    return(invisible(NULL))
  }

  # determine two subgroups (we assume exactly 2 for the main story)
  subs <- sort(unique(m[[sub_col]]))
  if (length(subs) < 2) stop2("Condition '%s' has <2 subgroups in %s.", cond, sub_col)
  if (length(subs) > 2) {
    log_msg("Condition '%s' has %d subgroups. Using first two (lexicographic): %s",
            cond, length(subs), paste(subs, collapse=", "))
    subs <- subs[1:2]
  }
  s1 <- subs[1]; s2 <- subs[2]

  # sample IDs
  treat_all <- m[get(grp_col) == "treatment", get(sid_col)]
  ctrl_all  <- m[get(grp_col) == "control",   get(sid_col)]

  m1 <- m[get(sub_col) == s1]
  m2 <- m[get(sub_col) == s2]

  treat_1 <- m1[get(grp_col) == "treatment", get(sid_col)]
  ctrl_1  <- m1[get(grp_col) == "control",   get(sid_col)]
  treat_2 <- m2[get(grp_col) == "treatment", get(sid_col)]
  ctrl_2  <- m2[get(grp_col) == "control",   get(sid_col)]

  log_msg("COND=%s | %s: treat=%d ctrl=%d | %s: treat=%d ctrl=%d | pooled: treat=%d ctrl=%d",
          cond, s1, length(treat_1), length(ctrl_1), s2, length(treat_2), length(ctrl_2),
          length(treat_all), length(ctrl_all))

  # effects
  d1 <- effect_delta(X, treat_1, ctrl_1)
  d2 <- effect_delta(X, treat_2, ctrl_2)
  d_pool <- effect_delta(X, treat_all, ctrl_all)

  # stratified weighted by subgroup sizes (treat+ctrl)
  n1 <- length(treat_1) + length(ctrl_1)
  n2 <- length(treat_2) + length(ctrl_2)
  w1 <- n1 / (n1 + n2)
  w2 <- n2 / (n1 + n2)
  d_strat <- w1 * d1 + w2 * d2

  # variance proxy across all samples in this condition
  cond_samples <- m[[sid_col]]
  var_all <- apply(X[, cond_samples, drop=FALSE], 1, var)

  # instability metrics
  eps <- 1e-6
  sgn1 <- vapply(d1, sgn, integer(1), eps=eps)
  sgn2 <- vapply(d2, sgn, integer(1), eps=eps)

  sign_flip <- (sgn1 != 0L) & (sgn2 != 0L) & (sgn1 != sgn2)
  hetero_gap <- abs(d1 - d2)
  discordance <- abs(d_pool - d_strat)

  # simple combined score (transparent)
  instability_score <- hetero_gap + discordance

  out <- data.table(
    feature = rownames(X),
    condition = cond,
    subgroup_1 = s1,
    subgroup_2 = s2,
    delta_1 = d1,
    delta_2 = d2,
    delta_pool = d_pool,
    delta_strat = d_strat,
    sign_flip = sign_flip,
    hetero_gap = hetero_gap,
    discordance = discordance,
    instability_score = instability_score,
    var_proxy = var_all
  )

  # rank tables
  out_rank <- out[order(-instability_score)]
  out_top <- out_rank[1:min(topn, .N)]

  # save
  out_path <- file.path(tables_dir, sprintf("%s_pool_vs_strat_all.tsv", cond))
  top_path <- file.path(tables_dir, sprintf("%s_pool_vs_strat_top%d.tsv", cond, topn))
  fwrite(out, out_path, sep="\t")
  fwrite(out_top, top_path, sep="\t")

  # summary
  sum_path <- file.path(tables_dir, sprintf("%s_pool_vs_strat_summary.tsv", cond))
  summary <- out[, .(
    n_features = .N,
    n_sign_flip = sum(sign_flip, na.rm=TRUE),
    frac_sign_flip = mean(sign_flip, na.rm=TRUE),
    median_hetero_gap = median(hetero_gap, na.rm=TRUE),
    median_discordance = median(discordance, na.rm=TRUE),
    median_instability = median(instability_score, na.rm=TRUE),
    cor_instability_vs_var = suppressWarnings(cor(instability_score, var_proxy, use="pairwise.complete.obs"))
  )]
  fwrite(summary, sum_path, sep="\t")

  log_msg("Saved: %s", out_path)
  log_msg("Saved: %s", top_path)
  log_msg("Saved: %s", sum_path)
}

# run for UV + DES (your current main story)
run_condition("uv")
run_condition("des")

log_msg("NM2 pool-vs-strat done.")
