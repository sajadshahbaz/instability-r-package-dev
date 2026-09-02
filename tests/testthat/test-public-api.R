test_that("public API has exactly the locked formal arguments", {
  expect_true(is.function(compute_uv_des_instability))
  expect_identical(names(formals(compute_uv_des_instability)), c("vst", "metadata"))
})

test_that("public API returns the locked ordinary eleven-table list", {
  vst <- data.frame(
    feature = c("gene_A", "gene_B", "gene_C"),
    uv_a_control = c(1, 2, 3),
    uv_a_treatment = c(2, 1, 5),
    uv_b_control = c(4, 3, 2),
    uv_b_treatment = c(2, 6, 3),
    des_a_control = c(2, 3, 4),
    des_a_treatment = c(4, 2, 7),
    des_b_control = c(5, 4, 3),
    des_b_treatment = c(3, 7, 4),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample_id = names(vst)[-1],
    condition = rep(c("uv", "des"), each = 4),
    baseline_block = rep(c("block_a", "block_a", "block_b", "block_b"), 2),
    group_type = rep(c("control", "treatment", "control", "treatment"), 2)
  )
  expected_names <- names(instability:::.authoritative_output_files())
  working_directory <- getwd()
  before <- list.files(tempdir(), pattern = "^instability-", full.names = TRUE)

  result <- compute_uv_des_instability(vst, metadata)

  after <- list.files(tempdir(), pattern = "^instability-", full.names = TRUE)
  expect_type(result, "list")
  expect_identical(class(result), "list")
  expect_length(result, 11L)
  expect_identical(names(result), expected_names)
  expect_true(all(vapply(result, data.table::is.data.table, logical(1))))
  expect_identical(result$uv_pool_vs_strat_all$feature, vst$feature)
  expect_identical(result$des_pool_vs_strat_all$feature, vst$feature)
  expect_identical(result$core_instability_signal_table$condition, c("des", "uv"))
  expect_identical(getwd(), working_directory)
  expect_identical(after, before)
})
