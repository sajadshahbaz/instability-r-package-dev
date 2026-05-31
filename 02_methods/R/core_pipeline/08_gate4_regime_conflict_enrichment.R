#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

in_dir  <- file.path(nm2_root, "04_results", "tables")
out_dir <- file.path(nm2_root, "04_results", "gate4_regime_conflict_enrichment")

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

clean_numeric <- function(x) {
  as.numeric(gsub("[^0-9eE+.-]", "", as.character(x)))
}

summarize_group <- function(dt, label) {
  data.table(
    group = label,
    n = nrow(dt),
    sign_flip_rate = mean(dt$sign_flip_logical, na.rm = TRUE),
    median_instability = median(dt$instability_score, na.rm = TRUE),
    median_hetero_gap = median(dt$hetero_gap, na.rm = TRUE),
    median_discordance = median(dt$discordance, na.rm = TRUE)
  )
}

run_condition <- function(cond) {
  cat("Processing:", cond, "\n")

  path <- file.path(in_dir, paste0(cond, "_pool_vs_strat_all.tsv"))
  dt <- fread(path)

  num_cols <- intersect(c("instability_score", "hetero_gap", "discordance", "var_proxy", "delta_pool"), names(dt))
  for (cc in num_cols) {
    dt[[cc]] <- clean_numeric(dt[[cc]])
  }

  if ("sign_flip" %in% names(dt)) {
    dt[, sign_flip_logical := as.character(sign_flip) %in% c("TRUE", "True", "true", "1")]
  } else if ("discordance" %in% names(dt)) {
    dt[, sign_flip_logical := discordance > 0]
  } else {
    dt[, sign_flip_logical := FALSE]
  }

  cut90 <- quantile(dt$instability_score, 0.90, na.rm = TRUE)
  cut95 <- quantile(dt$instability_score, 0.95, na.rm = TRUE)

  bg    <- summarize_group(dt, "background_all")
  top10 <- summarize_group(dt[instability_score >= cut90], "top10_instability")
  top5  <- summarize_group(dt[instability_score >= cut95], "top5_instability")

  summary_tab <- rbindlist(list(bg, top10, top5))
  summary_tab[, sign_flip_rate := round(sign_flip_rate, 4)]
  summary_tab[, median_instability := round(median_instability, 4)]
  summary_tab[, median_hetero_gap := round(median_hetero_gap, 4)]
  summary_tab[, median_discordance := round(median_discordance, 4)]
  summary_tab[, condition := cond]
  setcolorder(summary_tab, c("condition", "group", "n",
                             "sign_flip_rate", "median_instability",
                             "median_hetero_gap", "median_discordance"))

  fwrite(summary_tab,
         file.path(tab_dir, paste0(cond, "_regime_conflict_summary.tsv")),
         sep = "\t")

  sf_plot_dt <- copy(summary_tab)
  sf_plot_dt[, group := factor(group,
                               levels = c("background_all", "top10_instability", "top5_instability"))]

  p1 <- ggplot(sf_plot_dt, aes(x = group, y = sign_flip_rate)) +
    geom_col() +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(toupper(cond), ": sign-flip enrichment in high-instability genes"),
      x = NULL,
      y = "sign-flip rate"
    )

  ggsave(file.path(fig_dir, paste0(cond, "_signflip_enrichment.png")),
         p1, width = 7, height = 5, dpi = 300)

  dist_dt <- rbindlist(list(
    dt[, .(group = "background_all", hetero_gap)],
    dt[instability_score >= cut90, .(group = "top10_instability", hetero_gap)],
    dt[instability_score >= cut95, .(group = "top5_instability", hetero_gap)]
  ))

  dist_dt[, group := factor(group,
                            levels = c("background_all", "top10_instability", "top5_instability"))]

  p2 <- ggplot(dist_dt, aes(x = group, y = hetero_gap)) +
    geom_boxplot(outlier.alpha = 0.2) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(toupper(cond), ": hetero_gap distribution by instability tier"),
      x = NULL,
      y = "hetero_gap"
    )

  ggsave(file.path(fig_dir, paste0(cond, "_heterogap_by_instability.png")),
         p2, width = 7, height = 5, dpi = 300)

  id_col <- if ("feature" %in% names(dt)) "feature" else if ("gene" %in% names(dt)) "gene" else stop("No feature/gene column found")

  keep_cols <- intersect(
    c(id_col, "condition", "subgroup_1", "subgroup_2",
      "delta_1", "delta_2", "delta_pool", "delta_strat",
      "sign_flip", "hetero_gap", "discordance", "instability_score", "var_proxy", "n_timepoints"),
    names(dt)
  )

  top_features <- dt[instability_score >= cut95][order(-instability_score, -hetero_gap)][
    1:min(.N, 50), ..keep_cols
  ]

  fwrite(top_features,
         file.path(tab_dir, paste0(cond, "_top5_instability_features.tsv")),
         sep = "\t")

  return(summary_tab)
}

all_summary <- rbindlist(lapply(c("uv", "des", "gam"), run_condition))

fwrite(all_summary,
       file.path(tab_dir, "gate4_regime_conflict_summary_all.tsv"),
       sep = "\t")

cat("\nGate 4 completed.\n")
