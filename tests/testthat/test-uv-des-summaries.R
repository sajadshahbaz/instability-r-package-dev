test_that("safe Spearman uses finite pairs and requires three", {
  expect_identical(
    instability:::.safe_spearman(c(1, 2, 3, Inf, NA), c(3, 2, 1, 4, 5)),
    -1
  )
  expect_true(is.na(instability:::.safe_spearman(c(1, 2, Inf), c(1, 2, 3))))
})

test_that("Jaccard preserves unique-set and empty-union behavior", {
  expect_identical(instability:::.jaccard(c("a", "a", "b"), c("b", "c")), 1 / 3)
  expect_true(is.na(instability:::.jaccard(character(), character())))
})

test_that("tail tables use locked K values and preserve missingness", {
  dt <- data.table::data.table(
    feature = c("a", "b", "c"),
    hetero_gap = c(3, 2, 1),
    discordance = c(1, 2, 3),
    instability_score = c(4, 4, 4),
    delta_pool = c(1, NA_real_, 3),
    delta_strat = c(0, 0, 0),
    delta_1 = c(1, 2, 3),
    delta_2 = c(0, 0, 0),
    sign_flip = c(TRUE, FALSE, TRUE)
  )

  tail <- instability:::.make_tail_table(dt, "uv")

  expect_identical(tail$K, c(50, 100, 200, 500, 1000, 2000))
  expect_identical(tail$condition, rep("uv", 6))
  expect_true(all(is.na(tail$median_abs_pool_minus_strat_in_top_comb)))
  expect_identical(tail$j_gap_vs_dis, rep(1, 6))
})

test_that("rank vectors retain finite scores, ordering, and ties", {
  dt <- data.table::data.table(
    feature = c("a", "b", "c", "d"),
    score = c(2, NA_real_, 2, Inf)
  )

  rank <- instability:::.rank_vector(dt, "score")

  expect_identical(rank$feature, c("a", "c"))
  expect_identical(rank$rank, 1:2)
})

test_that("rank displacement preserves merge order, signs, and x", {
  dt <- data.table::data.table(
    feature = c("b", "a", "c"),
    hetero_gap = c(3, 2, 1),
    discordance = c(1, 3, 2)
  )

  displacement <- instability:::.make_rank_displacement(dt, "des")

  expect_identical(displacement$feature, c("a", "b", "c"))
  expect_identical(displacement$rank_strat, c(2L, 1L, 3L))
  expect_identical(displacement$rank_pool, c(1L, 3L, 2L))
  expect_identical(displacement$abs_rank_diff, c(1L, 2L, 1L))
  expect_identical(displacement$signed_rank_diff, c(-1L, 2L, -1L))
  expect_identical(displacement$condition, rep("des", 3))
  expect_identical(displacement$x, log10(1 + c(1L, 2L, 1L)))
  expect_identical(
    names(displacement),
    c("feature", "rank_strat", "rank_pool", "abs_rank_diff", "signed_rank_diff", "condition", "x")
  )
})

test_that("file stages preserve condition ordering and headline extraction", {
  directory <- tempfile("instability-summary-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  make_condition <- function(condition) {
    data.table::data.table(
      feature = sprintf("g%03d", seq_len(60)),
      condition = condition,
      subgroup_1 = "a",
      subgroup_2 = "b",
      delta_1 = seq_len(60),
      delta_2 = -seq_len(60),
      delta_pool = seq_len(60) / 10,
      delta_strat = seq_len(60) / 20,
      sign_flip = rep(c(TRUE, FALSE), 30),
      hetero_gap = seq_len(60),
      discordance = rev(seq_len(60)),
      instability_score = seq_len(60) + rev(seq_len(60)),
      var_proxy = seq_len(60)
    )
  }

  data.table::fwrite(make_condition("uv"), file.path(directory, "uv_pool_vs_strat_all.tsv"), sep = "\t")
  data.table::fwrite(make_condition("des"), file.path(directory, "des_pool_vs_strat_all.tsv"), sep = "\t")

  metrics <- instability:::.make_pool_strat_metrics(directory)
  tails <- instability:::.make_tail_outputs(directory)
  core <- instability:::.make_core_table(directory)

  expect_identical(metrics$condition, c("des", "uv"))
  expect_identical(tails$tail_headline_summary$condition, c("des", "uv"))
  expect_identical(
    names(tails$tail_headline_summary),
    c("condition", "j_gap_vs_dis_K200", "j_gap_vs_dis_K1000", "flip_top200", "flip_top1000")
  )
  expect_identical(core$condition, c("des", "uv"))
  expect_identical(
    names(core),
    c(
      "condition", "n_features", "spearman_instability_vs_effect_divergence",
      "spearman_pool_vs_strat", "median_abs_pool_minus_strat",
      "q95_abs_pool_minus_strat", "sign_flip_rate", "j_gap_vs_dis_K200",
      "j_gap_vs_dis_K1000", "flip_top200", "flip_top1000"
    )
  )
})
