test_that("metadata normalization and sample intersection follow the source", {
  vst <- data.table::data.table(feature = "g1", s2 = 20, s1 = 10)
  metadata <- data.table::data.table(
    sample_id = c(" s1 ", "s2", "outside"),
    condition = c(" UV ", "uv", "DES"),
    baseline_block = c(" Long ", "SHORT", "fast"),
    group_type = c(" Control ", "TREATMENT", "control")
  )

  prepared <- instability:::.prepare_uv_des_inputs(vst, metadata)

  expect_identical(colnames(prepared$X), c("s1", "s2"))
  expect_identical(unname(prepared$X[1, ]), c(10, 20))
  expect_identical(prepared$meta$sample_id, c("s1", "s2"))
  expect_identical(prepared$meta$condition, c("uv", "uv"))
  expect_identical(prepared$meta$baseline_block, c("long", "short"))
  expect_identical(prepared$meta$group_type, c("control", "treatment"))
})

test_that("effect deltas preserve rowMeans and empty-group behavior", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2, dimnames = list(c("g1", "g2"), c("t", "c")))

  expect_identical(instability:::.effect_delta(mat, "t", "c"), c(g1 = -3, g2 = -6))
  expect_true(all(is.na(instability:::.effect_delta(mat, character(), "c"))))
  expect_true(all(is.na(instability:::.effect_delta(mat, "t", character()))))

  mat[1, "t"] <- NA_real_
  expect_true(is.na(instability:::.effect_delta(mat, "t", "c")[[1]]))
})

test_that("epsilon sign behavior is scalar and inclusive at the boundary", {
  values <- c(NA_real_, -2e-6, -1e-6, 0, 1e-6, 2e-6)
  observed <- vapply(values, instability:::.sign_with_epsilon, integer(1), eps = 1e-6)
  expect_identical(observed, c(NA_integer_, -1L, 0L, 0L, 0L, 1L))
})

test_that("condition computation uses the first two lexicographic subgroups", {
  samples <- c("a_c", "a_t", "b_c", "b_t", "c_c", "c_t")
  X <- matrix(
    c(0, 2, 10, 6, 100, 130),
    nrow = 1,
    dimnames = list("g1", samples)
  )
  meta <- data.table::data.table(
    sample_id = samples,
    condition = "uv",
    baseline_block = c("a", "a", "b", "b", "c", "c"),
    group_type = rep(c("control", "treatment"), 3)
  )
  data.table::setkey(meta, sample_id)

  result <- instability:::.compute_uv_des_condition(X, meta, "uv")

  expect_identical(result$subgroup_1, "a")
  expect_identical(result$subgroup_2, "b")
  expect_identical(result$delta_1, 2)
  expect_identical(result$delta_2, -4)
  expect_identical(result$delta_pool, 28 / 3)
  expect_identical(result$delta_strat, -1)
  expect_identical(result$sign_flip, TRUE)
  expect_identical(result$hetero_gap, 6)
  expect_identical(result$discordance, 31 / 3)
  expect_identical(result$instability_score, 49 / 3)
})

test_that("stratified effects use treatment-plus-control subgroup weights", {
  samples <- c("a_c1", "a_t1", "b_c1", "b_c2", "b_t1", "b_t2")
  X <- matrix(
    c(0, 4, 10, 14, 12, 16),
    nrow = 1,
    dimnames = list("g1", samples)
  )
  meta <- data.table::data.table(
    sample_id = samples,
    condition = "des",
    baseline_block = c("a", "a", "b", "b", "b", "b"),
    group_type = c("control", "treatment", "control", "control", "treatment", "treatment")
  )
  data.table::setkey(meta, sample_id)

  result <- instability:::.compute_uv_des_condition(X, meta, "des")

  expect_identical(result$delta_1, 4)
  expect_identical(result$delta_2, 2)
  expect_identical(result$delta_strat, 8 / 3)
})

test_that("Top-150 ordering retains source tie order", {
  table <- data.table::data.table(
    feature = sprintf("g%03d", seq_len(151)),
    instability_score = c(10, 10, rev(seq_len(149)))
  )

  top <- instability:::.make_top150(table)

  expect_identical(nrow(top), 150L)
  expect_identical(top$feature[1:2], c("g003", "g004"))
  tie_rows <- top[instability_score == 10, feature]
  expect_identical(tie_rows, c("g141", "g001", "g002"))
})

test_that("condition summary preserves missing-value and correlation behavior", {
  table <- data.table::data.table(
    sign_flip = c(TRUE, FALSE, NA),
    hetero_gap = c(1, 3, NA_real_),
    discordance = c(2, 4, NA_real_),
    instability_score = c(3, 7, NA_real_),
    var_proxy = c(6, 14, 20)
  )

  summary <- instability:::.summarize_uv_des_condition(table)

  expect_identical(summary$n_features, 3L)
  expect_identical(summary$n_sign_flip, 1L)
  expect_identical(summary$frac_sign_flip, 0.5)
  expect_identical(summary$median_hetero_gap, 2)
  expect_identical(summary$median_discordance, 3)
  expect_identical(summary$median_instability, 5)
  expect_identical(summary$cor_instability_vs_var, 1)
})
