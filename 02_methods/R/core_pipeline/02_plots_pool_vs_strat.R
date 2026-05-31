#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- Sys.getenv("NM2_ROOT")
if (ROOT == "") {
  ROOT <- normalizePath(getwd(), mustWork = TRUE)
  Sys.setenv(NM2_ROOT = ROOT)
}

# Warm up graphics device to reduce fontconfig noise on some systems
try(suppressWarnings(grDevices::pdf(NULL)), silent = TRUE)
try(suppressWarnings(grDevices::dev.off()), silent = TRUE)

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
  cat(sprintf("[%s] ", ts()), sprintf(...), "\n", sep = "")
}

stop2 <- function(...) {
  cat(sprintf("[FATAL %s] ", ts()), sprintf(...), "\n", sep = "", file = stderr())
  quit(status = 2)
}

out_fig <- out_fig <- file.path(ROOT, "04_results", "qc_figures", "pool_vs_strat")
out_tab <- file.path(ROOT, "04_results", "tables")
out_log <- file.path(ROOT, "04_results", "logs")

dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log, recursive = TRUE, showWarnings = FALSE)

log10p1 <- function(x) log10(1 + pmax(0, x))

need_cols <- c(
  "feature",
  "hetero_gap",
  "discordance",
  "instability_score",
  "delta_pool",
  "delta_strat",
  "delta_1",
  "delta_2",
  "sign_flip"
)

check_cols <- function(dt, cond) {
  miss <- setdiff(need_cols, names(dt))
  if (length(miss) > 0) {
    stop2("[%s] Missing columns: %s", cond, paste(miss, collapse = ", "))
  }
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

plot_instability_vs_effect_divergence <- function(dt, cond) {
  dt[, effect_divergence := abs(delta_1 - delta_2)]
  dt[, x := log10p1(instability_score)]
  dt[, y := log10p1(effect_divergence)]

  p <- ggplot(dt, aes(x = x, y = y)) +
    geom_point(size = 0.6, alpha = 0.35) +
    labs(
      title = sprintf("%s: Instability vs subgroup effect divergence", toupper(cond)),
      x = "log10(1 + instability score)",
      y = "log10(1 + |delta_1 - delta_2|)"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(out_fig, sprintf("%s_instability_vs_effect_divergence.pdf", cond)),
    p,
    width = 7,
    height = 5.2,
    dpi = 300
  )
}

plot_pooled_vs_stratified <- function(dt, cond) {
  dt[, x := delta_strat]
  dt[, y := delta_pool]

  p <- ggplot(dt, aes(x = x, y = y)) +
    geom_point(size = 0.6, alpha = 0.35) +
    labs(
      title = sprintf("%s: Pooled vs stratified effect estimates", toupper(cond)),
      x = "delta_strat (stratified estimate)",
      y = "delta_pool (pooled estimate)"
    ) +
    theme_classic(base_size = 12)

  ggsave(
  file.path(out_fig, sprintf("%s_pooled_vs_stratified.pdf", cond)),
  p,
  width = 7,
  height = 5.2,
  dpi = 300,
  device = "pdf"
)
}

main <- function() {
  log_file <- file.path(out_log, "02_plots_pool_vs_strat_log.txt")
  sink(log_file, split = TRUE)
  on.exit(sink(), add = TRUE)

  log_msg("NM2 plotting start")
  log_msg("NM2_ROOT = %s", ROOT)
  log_msg("Input table directory = %s", out_tab)
  log_msg("Output figure directory = %s", out_fig)

  files <- list.files(
    out_tab,
    pattern = "_pool_vs_strat_all\\.tsv$",
    full.names = TRUE
  )
files <- files[grepl("(uv|des)_pool_vs_strat_all\\.tsv$", basename(files))]
  if (length(files) == 0) {
  stop2("No UV/DES *_pool_vs_strat_all.tsv files found in: %s", out_tab)
}

  metrics <- list()

  for (f in files) {
    cond <- sub("_pool_vs_strat_all\\.tsv$", "", basename(f))
    log_msg("Processing condition: %s", cond)

    dt <- fread(f)
    check_cols(dt, cond)

    dt[, effect_divergence := abs(delta_1 - delta_2)]

    m <- data.table(
      condition = cond,
      n_features = nrow(dt),
      spearman_instability_vs_effect_divergence =
        safe_spearman(dt$instability_score, dt$effect_divergence),
      spearman_pool_vs_strat =
        safe_spearman(dt$delta_pool, dt$delta_strat),
      median_abs_pool_minus_strat =
        median(abs(dt$delta_pool - dt$delta_strat), na.rm = TRUE),
      q95_abs_pool_minus_strat =
        as.numeric(quantile(abs(dt$delta_pool - dt$delta_strat), 0.95, na.rm = TRUE)),
      sign_flip_rate =
        mean(as.numeric(dt$sign_flip), na.rm = TRUE)
    )

    metrics[[length(metrics) + 1]] <- m

    plot_instability_vs_effect_divergence(dt, cond)
    plot_pooled_vs_stratified(dt, cond)
  }

  metrics_dt <- rbindlist(metrics)
  setorder(metrics_dt, condition)

  out_metrics <- file.path(out_tab, "pool_vs_strat_metrics.tsv")
  fwrite(metrics_dt, out_metrics, sep = "\t")

  log_msg("Saved metrics: %s", out_metrics)
  log_msg("NM2 plotting done")
}

main()
