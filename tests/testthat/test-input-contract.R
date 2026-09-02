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

test_that("public validation rejects structurally empty inputs actionably", {
  metadata <- data.frame(
    sample_id = "s1",
    condition = "uv",
    baseline_block = "a",
    group_type = "control"
  )

  expect_error(
    compute_uv_des_instability(data.frame(), metadata),
    "vst must contain a first feature-identifier column and at least one sample column",
    fixed = TRUE
  )
  expect_error(
    compute_uv_des_instability(data.frame(feature = character(), s1 = numeric()), metadata),
    "vst must contain at least one feature row",
    fixed = TRUE
  )
  expect_error(
    compute_uv_des_instability(data.frame(feature = "g1"), metadata),
    "vst must contain at least one named sample column",
    fixed = TRUE
  )
  expect_error(
    compute_uv_des_instability(
      data.frame(feature = "g1", s1 = 1),
      metadata[FALSE, ]
    ),
    "metadata must contain at least one sample row",
    fixed = TRUE
  )
})

test_that("public validation requires matched UV and DES conditions", {
  vst <- data.frame(
    feature = c("g1", "g2"),
    uv_a_c = c(1, 2),
    uv_a_t = c(2, 3),
    uv_b_c = c(3, 4),
    uv_b_t = c(4, 5),
    des_a_c = c(1, 3),
    des_a_t = c(3, 5),
    des_b_c = c(2, 4),
    des_b_t = c(4, 6),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample_id = names(vst)[-1],
    condition = rep(c("uv", "des"), each = 4),
    baseline_block = rep(c("a", "a", "b", "b"), 2),
    group_type = rep(c("control", "treatment", "control", "treatment"), 2)
  )

  expect_error(
    compute_uv_des_instability(vst, metadata[metadata$condition == "uv", ]),
    "metadata must contain matched samples for required condition(s): des",
    fixed = TRUE
  )
  expect_error(
    compute_uv_des_instability(vst, metadata[metadata$condition == "des", ]),
    "metadata must contain matched samples for required condition(s): uv",
    fixed = TRUE
  )
})
