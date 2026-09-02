test_that("public API has exactly the locked formal arguments", {
  expect_true(is.function(compute_uv_des_instability))
  expect_identical(names(formals(compute_uv_des_instability)), c("vst", "metadata"))
})

test_that("public API returns the locked ordinary eleven-table list", {
  samples <- c(
    "uv_a_c", "uv_a_t", "uv_b_c", "uv_b_t",
    "des_a_c", "des_a_t", "des_b_c", "des_b_t"
  )
  vst <- data.table::data.table(
    feature = c("g1", "g2", "g3"),
    uv_a_c = c(0, 2, 4),
    uv_a_t = c(2, 1, 8),
    uv_b_c = c(4, 3, 2),
    uv_b_t = c(1, 6, 3),
    des_a_c = c(1, 3, 5),
    des_a_t = c(3, 2, 9),
    des_b_c = c(5, 4, 3),
    des_b_t = c(2, 7, 4)
  )
  metadata <- data.table::data.table(
    sample_id = samples,
    condition = rep(c("uv", "des"), each = 4),
    baseline_block = rep(c("a", "a", "b", "b"), 2),
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
  expect_identical(getwd(), working_directory)
  expect_identical(after, before)
})
