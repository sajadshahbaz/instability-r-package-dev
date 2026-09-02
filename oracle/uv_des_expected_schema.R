uv_des_expected <- list(
  scientific_files = c(
    "uv_pool_vs_strat_all.tsv",
    "uv_pool_vs_strat_top150.tsv",
    "uv_pool_vs_strat_summary.tsv",
    "des_pool_vs_strat_all.tsv",
    "des_pool_vs_strat_top150.tsv",
    "des_pool_vs_strat_summary.tsv",
    "pool_vs_strat_metrics.tsv",
    "tail_jaccard_metrics.tsv",
    "rank_displacement_pool_vs_strat.tsv",
    "tail_headline_summary.tsv",
    "core_instability_signal_table.tsv"
  ),
  schemas = list(
    pool_all = c(
      "feature", "condition", "subgroup_1", "subgroup_2", "delta_1",
      "delta_2", "delta_pool", "delta_strat", "sign_flip", "hetero_gap",
      "discordance", "instability_score", "var_proxy"
    ),
    pool_summary = c(
      "n_features", "n_sign_flip", "frac_sign_flip", "median_hetero_gap",
      "median_discordance", "median_instability", "cor_instability_vs_var"
    ),
    metrics = c(
      "condition", "n_features", "spearman_instability_vs_effect_divergence",
      "spearman_pool_vs_strat", "median_abs_pool_minus_strat",
      "q95_abs_pool_minus_strat", "sign_flip_rate"
    ),
    tail = c(
      "condition", "K", "j_gap_vs_dis", "j_gap_vs_comb", "j_dis_vs_comb",
      "flip_rate_in_top_comb", "median_abs_pool_minus_strat_in_top_comb",
      "median_hetero_gap_in_top_comb", "median_discordance_in_top_comb"
    ),
    rank = c(
      "feature", "rank_strat", "rank_pool", "abs_rank_diff",
      "signed_rank_diff", "condition", "x"
    ),
    headline = c(
      "condition", "j_gap_vs_dis_K200", "j_gap_vs_dis_K1000",
      "flip_top200", "flip_top1000"
    ),
    core = c(
      "condition", "n_features", "spearman_instability_vs_effect_divergence",
      "spearman_pool_vs_strat", "median_abs_pool_minus_strat",
      "q95_abs_pool_minus_strat", "sign_flip_rate", "j_gap_vs_dis_K200",
      "j_gap_vs_dis_K1000", "flip_top200", "flip_top1000"
    )
  ),
  dimensions = list(
    uv_pool_vs_strat_all.tsv = c(15187L, 13L),
    uv_pool_vs_strat_top150.tsv = c(150L, 13L),
    uv_pool_vs_strat_summary.tsv = c(1L, 7L),
    des_pool_vs_strat_all.tsv = c(15187L, 13L),
    des_pool_vs_strat_top150.tsv = c(150L, 13L),
    des_pool_vs_strat_summary.tsv = c(1L, 7L),
    pool_vs_strat_metrics.tsv = c(2L, 7L),
    tail_jaccard_metrics.tsv = c(12L, 9L),
    rank_displacement_pool_vs_strat.tsv = c(30374L, 7L),
    tail_headline_summary.tsv = c(2L, 5L),
    core_instability_signal_table.tsv = c(2L, 11L)
  ),
  k_values = c(50L, 100L, 200L, 500L, 1000L, 2000L),
  condition_order = c("des", "uv")
)

uv_des_schema_for <- function(filename) {
  if (grepl("_pool_vs_strat_(all|top150)\\.tsv$", filename)) {
    return(uv_des_expected$schemas$pool_all)
  }
  if (grepl("_pool_vs_strat_summary\\.tsv$", filename)) {
    return(uv_des_expected$schemas$pool_summary)
  }
  switch(filename,
    pool_vs_strat_metrics.tsv = uv_des_expected$schemas$metrics,
    tail_jaccard_metrics.tsv = uv_des_expected$schemas$tail,
    rank_displacement_pool_vs_strat.tsv = uv_des_expected$schemas$rank,
    tail_headline_summary.tsv = uv_des_expected$schemas$headline,
    core_instability_signal_table.tsv = uv_des_expected$schemas$core,
    stop("No expected schema for: ", filename, call. = FALSE)
  )
}
