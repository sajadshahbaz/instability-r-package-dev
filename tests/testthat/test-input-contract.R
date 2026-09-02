test_that("Stage-2A input contract accepts only tabular inputs", {
  metadata <- data.frame(
    sample_id = "s1",
    condition = "uv",
    baseline_block = "a",
    group_type = "control"
  )
  vst <- data.frame(feature = "g1", s1 = 1, check.names = FALSE)

  expect_error(
    instability:::.prepare_uv_des_inputs(as.matrix(vst), metadata),
    "vst must be a data.frame or data.table",
    fixed = TRUE
  )
  expect_error(
    instability:::.prepare_uv_des_inputs(vst, as.matrix(metadata)),
    "metadata must be a data.frame or data.table",
    fixed = TRUE
  )
})

test_that("Stage-2A input contract requires canonical metadata columns", {
  vst <- data.frame(feature = "g1", s1 = 1, check.names = FALSE)
  metadata <- data.frame(
    sample_id = "s1",
    condition = "uv",
    baseline_block = "a"
  )

  expect_error(
    instability:::.prepare_uv_des_inputs(vst, metadata),
    "metadata missing required columns: group_type",
    fixed = TRUE
  )
})
