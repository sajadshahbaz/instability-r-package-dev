.safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

.make_pool_strat_metrics <- function(table_directory) {
  files <- list.files(
    table_directory,
    pattern = "_pool_vs_strat_all\\.tsv$",
    full.names = TRUE
  )
  files <- files[grepl("(uv|des)_pool_vs_strat_all\\.tsv$", basename(files))]
  if (length(files) == 0) {
    stop(sprintf("No UV/DES *_pool_vs_strat_all.tsv files found in: %s", table_directory), call. = FALSE)
  }

  metrics <- list()

  for (f in files) {
    cond <- sub("_pool_vs_strat_all\\.tsv$", "", basename(f))
    dt <- data.table::fread(f)

    need_cols <- c(
      "feature", "hetero_gap", "discordance", "instability_score",
      "delta_pool", "delta_strat", "delta_1", "delta_2", "sign_flip"
    )
    miss <- setdiff(need_cols, names(dt))
    if (length(miss) > 0) {
      stop(sprintf("[%s] Missing columns: %s", cond, paste(miss, collapse = ", ")), call. = FALSE)
    }

    dt[, effect_divergence := abs(delta_1 - delta_2)]

    m <- data.table::data.table(
      condition = cond,
      n_features = nrow(dt),
      spearman_instability_vs_effect_divergence =
        .safe_spearman(dt$instability_score, dt$effect_divergence),
      spearman_pool_vs_strat =
        .safe_spearman(dt$delta_pool, dt$delta_strat),
      median_abs_pool_minus_strat =
        median(abs(dt$delta_pool - dt$delta_strat), na.rm = TRUE),
      q95_abs_pool_minus_strat =
        as.numeric(quantile(abs(dt$delta_pool - dt$delta_strat), 0.95, na.rm = TRUE)),
      sign_flip_rate =
        mean(as.numeric(dt$sign_flip), na.rm = TRUE)
    )

    metrics[[length(metrics) + 1]] <- m
  }

  metrics_dt <- data.table::rbindlist(metrics)
  data.table::setorder(metrics_dt, condition)

  data.table::fwrite(
    metrics_dt,
    file.path(table_directory, "pool_vs_strat_metrics.tsv"),
    sep = "\t"
  )

  metrics_dt
}

.jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

.make_tail_table <- function(dt, condition) {
  cond <- condition
  ks <- c(50, 100, 200, 500, 1000, 2000)
  out <- list()
  for (K in ks) {
    top_gap <- dt[order(-hetero_gap)][1:min(K, .N), feature]
    top_dis <- dt[order(-discordance)][1:min(K, .N), feature]
    top_comb <- dt[order(-instability_score)][1:min(K, .N), feature]

    out[[length(out) + 1]] <- data.table::data.table(
      condition = cond,
      K = K,
      j_gap_vs_dis = .jaccard(top_gap, top_dis),
      j_gap_vs_comb = .jaccard(top_gap, top_comb),
      j_dis_vs_comb = .jaccard(top_dis, top_comb),
      flip_rate_in_top_comb = dt[feature %in% top_comb, mean(as.numeric(sign_flip))],
      median_abs_pool_minus_strat_in_top_comb =
        dt[feature %in% top_comb, median(abs(delta_pool - delta_strat))],
      median_hetero_gap_in_top_comb = dt[feature %in% top_comb, median(hetero_gap)],
      median_discordance_in_top_comb = dt[feature %in% top_comb, median(discordance)]
    )
  }
  data.table::rbindlist(out)
}

.rank_vector <- function(dt, score_column) {
  if (!("feature" %in% names(dt))) {
    stop("rank_vec: missing column feature", call. = FALSE)
  }
  if (!(score_column %in% names(dt))) {
    stop(sprintf("rank_vec: missing score column: %s", score_column), call. = FALSE)
  }
  x <- dt[is.finite(get(score_column)), .(feature, score = get(score_column))]
  data.table::setorderv(x, cols = "score", order = -1L)
  x[, rank := seq_len(.N)]
  x[, .(feature, rank)]
}

.make_rank_displacement <- function(dt, condition) {
  cond <- condition
  r_strat <- .rank_vector(dt, "hetero_gap")
  r_pool <- .rank_vector(dt, "discordance")
  m <- merge(r_strat, r_pool, by = "feature", suffixes = c("_strat", "_pool"))
  m[, abs_rank_diff := abs(rank_strat - rank_pool)]
  m[, signed_rank_diff := (rank_pool - rank_strat)]
  m[, condition := cond]
  m[, x := log10(1 + abs_rank_diff)]
  m
}

.make_tail_outputs <- function(table_directory) {
  files <- list.files(
    table_directory,
    pattern = "_pool_vs_strat_all\\.tsv$",
    full.names = TRUE
  )
  files <- files[grepl("(uv|des)_pool_vs_strat_all\\.tsv$", basename(files))]
  if (length(files) == 0) {
    stop(sprintf("No UV/DES *_pool_vs_strat_all.tsv files found in: %s", table_directory), call. = FALSE)
  }

  tabs <- list()
  rds <- list()

  for (f in files) {
    cond <- sub("_pool_vs_strat_all\\.tsv$", "", basename(f))
    dt <- data.table::fread(f)

    need <- c(
      "feature", "hetero_gap", "discordance", "instability_score",
      "delta_pool", "delta_strat", "delta_1", "delta_2", "sign_flip"
    )
    miss <- setdiff(need, names(dt))
    if (length(miss) > 0) {
      stop(sprintf("[%s] Missing columns: %s", cond, paste(miss, collapse = ", ")), call. = FALSE)
    }

    tabs[[length(tabs) + 1]] <- .make_tail_table(dt, cond)
    rds[[length(rds) + 1]] <- .make_rank_displacement(dt, cond)
  }

  tab_all <- data.table::rbindlist(tabs)
  rd_all <- data.table::rbindlist(rds)

  data.table::fwrite(
    tab_all,
    file.path(table_directory, "tail_jaccard_metrics.tsv"),
    sep = "\t"
  )
  data.table::fwrite(
    rd_all,
    file.path(table_directory, "rank_displacement_pool_vs_strat.tsv"),
    sep = "\t"
  )

  headline <- tab_all[, .(
    j_gap_vs_dis_K200 = j_gap_vs_dis[K == 200][1],
    j_gap_vs_dis_K1000 = j_gap_vs_dis[K == 1000][1],
    flip_top200 = flip_rate_in_top_comb[K == 200][1],
    flip_top1000 = flip_rate_in_top_comb[K == 1000][1]
  ), by = condition]

  data.table::fwrite(
    headline,
    file.path(table_directory, "tail_headline_summary.tsv"),
    sep = "\t"
  )

  list(
    tail_jaccard_metrics = tab_all,
    rank_displacement_pool_vs_strat = rd_all,
    tail_headline_summary = headline
  )
}

.make_core_table <- function(table_directory) {
  in_metrics <- file.path(table_directory, "pool_vs_strat_metrics.tsv")
  in_headline <- file.path(table_directory, "tail_headline_summary.tsv")

  if (!file.exists(in_metrics)) {
    stop(sprintf("Missing: %s", in_metrics), call. = FALSE)
  }
  if (!file.exists(in_headline)) {
    stop(sprintf("Missing: %s", in_headline), call. = FALSE)
  }

  m <- data.table::fread(in_metrics)
  h <- data.table::fread(in_headline)

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
    stop(sprintf("Core table missing cols: %s", paste(miss, collapse = ", ")), call. = FALSE)
  }

  x <- x[, ..keep]
  data.table::setorder(x, condition)

  data.table::fwrite(
    x,
    file.path(table_directory, "core_instability_signal_table.tsv"),
    sep = "\t"
  )

  x
}
