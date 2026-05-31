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
try(suppressWarnings(grDevices::pdf(NULL)), silent=TRUE)
try(suppressWarnings(grDevices::dev.off()), silent=TRUE)

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) cat(sprintf("[%s] ", ts()), sprintf(...), "\n", sep="")
stop2 <- function(...) { cat(sprintf("[FATAL %s] ", ts()), sprintf(...), "\n", sep="", file=stderr()); quit(status=2) }

out_fig <- file.path(ROOT, "04_results", "figures")
out_tab <- file.path(ROOT, "04_results", "tables")
out_log <- file.path(ROOT, "04_results", "logs")
dir.create(out_fig, showWarnings=FALSE, recursive=TRUE)
dir.create(out_tab, showWarnings=FALSE, recursive=TRUE)
dir.create(out_log, showWarnings=FALSE, recursive=TRUE)

need <- c("feature","hetero_gap","discordance","instability_score","delta_pool","delta_strat","delta_1","delta_2","sign_flip")

check_cols <- function(dt, cond) {
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) stop2("[%s] Missing columns: %s", cond, paste(miss, collapse=", "))
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  inter <- length(intersect(a,b))
  uni <- length(union(a,b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

tail_table <- function(dt, cond, ks=c(50,100,200,500,1000,2000)) {
  out <- list()
  for (K in ks) {
    top_gap  <- dt[order(-hetero_gap)][1:min(K, .N), feature]
    top_dis  <- dt[order(-discordance)][1:min(K, .N), feature]
    top_comb <- dt[order(-instability_score)][1:min(K, .N), feature]

    out[[length(out)+1]] <- data.table(
      condition = cond,
      K = K,
      j_gap_vs_dis   = jaccard(top_gap, top_dis),
      j_gap_vs_comb  = jaccard(top_gap, top_comb),
      j_dis_vs_comb  = jaccard(top_dis, top_comb),
      flip_rate_in_top_comb = dt[feature %in% top_comb, mean(as.numeric(sign_flip))],
      median_abs_pool_minus_strat_in_top_comb =
        dt[feature %in% top_comb, median(abs(delta_pool - delta_strat))],
      median_hetero_gap_in_top_comb = dt[feature %in% top_comb, median(hetero_gap)],
      median_discordance_in_top_comb = dt[feature %in% top_comb, median(discordance)]
    )
  }
  rbindlist(out)
}

rank_vec <- function(dt, score_col) {
  if (!("feature" %in% names(dt))) stop2("rank_vec: missing column feature")
  if (!(score_col %in% names(dt))) stop2("rank_vec: missing score column: %s", score_col)
  x <- dt[is.finite(get(score_col)), .(feature, score = get(score_col))]
  setorderv(x, cols="score", order=-1L)  # descending
  x[, rank := seq_len(.N)]
  x[, .(feature, rank)]
}

rank_displacement <- function(dt, cond) {
  r_strat <- rank_vec(dt, "hetero_gap")
  r_pool  <- rank_vec(dt, "discordance")
  m <- merge(r_strat, r_pool, by="feature", suffixes=c("_strat","_pool"))
  m[, abs_rank_diff := abs(rank_strat - rank_pool)]
  m[, signed_rank_diff := (rank_pool - rank_strat)]
  m[, condition := cond]
  m
}

plot_jaccard_curve <- function(tab, cond) {
  pdt <- melt(
    tab,
    id.vars=c("condition","K"),
    measure.vars=c("j_gap_vs_dis","j_gap_vs_comb","j_dis_vs_comb"),
    variable.name="comparison",
    value.name="jaccard"
  )

  p <- ggplot(pdt, aes(x=K, y=jaccard, color=comparison, group=comparison)) +
    geom_line() + geom_point(size=1) +
    ylim(0, 1) +
    labs(
      title = sprintf("%s: Tail-set Jaccard vs K", toupper(cond)),
      x = "Top-K unstable features",
      y = "Jaccard overlap"
    ) +
    theme_classic(base_size=12)

  ggsave(file.path(out_fig, sprintf("%s_tail_jaccard_curve.png", cond)),
         p, width=7, height=5.2, dpi=300)
}

plot_rank_displacement <- function(rd, cond) {
  # allow zeros safely
  rd[, x := log10(1 + abs_rank_diff)]

  p <- ggplot(rd, aes(x=x)) +
    geom_histogram(bins=80) +
    labs(
      title = sprintf("%s: Rank displacement (pool vs strat)", toupper(cond)),
      x = "log10(1 + |rank(discordance) - rank(hetero_gap)|)",
      y = "Feature count"
    ) +
    theme_classic(base_size=12)

  ggsave(file.path(out_fig, sprintf("%s_rank_displacement_hist.png", cond)),
         p, width=7, height=5.2, dpi=300)
}

main <- function() {
  log_file <- file.path(out_log, "03_tail_divergence_log.txt")
  sink(log_file, split=TRUE); on.exit(sink(), add=TRUE)

  log_msg("NM2 tail divergence start")

  files <- list.files(out_tab, pattern="_pool_vs_strat_all\\.tsv$", full.names=TRUE)
  files <- files[grepl("(uv|des)_pool_vs_strat_all\\.tsv$", basename(files))]
  if (length(files) == 0) {
  stop2("No UV/DES *_pool_vs_strat_all.tsv files found in: %s", out_tab)
}

  ks <- c(50,100,200,500,1000,2000)

  tabs <- list()
  rds  <- list()

  for (f in files) {
    cond <- sub("_pool_vs_strat_all\\.tsv$", "", basename(f))
    dt <- fread(f)
    check_cols(dt, cond)

    tabs[[length(tabs)+1]] <- tail_table(dt, cond, ks)
    rds[[length(rds)+1]]   <- rank_displacement(dt, cond)

    # plots per condition
    plot_jaccard_curve(tabs[[length(tabs)]], cond)
    plot_rank_displacement(rds[[length(rds)]], cond)
  }

  tab_all <- rbindlist(tabs)
  rd_all  <- rbindlist(rds)

  out_tab1 <- file.path(out_tab, "tail_jaccard_metrics.tsv")
  fwrite(tab_all, out_tab1, sep="\t")
  log_msg("Saved: %s", out_tab1)

  out_tab2 <- file.path(out_tab, "rank_displacement_pool_vs_strat.tsv")
  fwrite(rd_all, out_tab2, sep="\t")
  log_msg("Saved: %s", out_tab2)

  # Headline summary (K=200 and K=1000) per condition
  headline <- tab_all[, .(
    j_gap_vs_dis_K200   = j_gap_vs_dis[K==200][1],
    j_gap_vs_dis_K1000  = j_gap_vs_dis[K==1000][1],
    flip_top200         = flip_rate_in_top_comb[K==200][1],
    flip_top1000        = flip_rate_in_top_comb[K==1000][1]
  ), by=condition]

  out_head <- file.path(out_tab, "tail_headline_summary.tsv")
  fwrite(headline, out_head, sep="\t")
  log_msg("Saved: %s", out_head)

  log_msg("NM2 tail divergence done")
}

main()
