test_that("authoritative output mapping has the locked names and order", {
  files <- instability:::.authoritative_output_files()

  expect_identical(
    names(files),
    c(
      "uv_pool_vs_strat_all",
      "uv_pool_vs_strat_top150",
      "uv_pool_vs_strat_summary",
      "des_pool_vs_strat_all",
      "des_pool_vs_strat_top150",
      "des_pool_vs_strat_summary",
      "pool_vs_strat_metrics",
      "tail_jaccard_metrics",
      "rank_displacement_pool_vs_strat",
      "tail_headline_summary",
      "core_instability_signal_table"
    )
  )
  expect_identical(unname(files), paste0(names(files), ".tsv"))
})

test_that("condition output writer uses authoritative filenames", {
  workspace <- instability:::.new_private_workspace()
  on.exit(instability:::.cleanup_private_workspace(workspace), add = TRUE)
  table <- data.table::data.table(value = 1)

  instability:::.write_condition_outputs(
    workspace,
    "uv",
    table,
    table,
    table
  )

  expect_true(file.exists(file.path(workspace, "uv_pool_vs_strat_all.tsv")))
  expect_true(file.exists(file.path(workspace, "uv_pool_vs_strat_top150.tsv")))
  expect_true(file.exists(file.path(workspace, "uv_pool_vs_strat_summary.tsv")))
})

test_that("final reread returns eleven data.tables in locked order", {
  workspace <- instability:::.new_private_workspace()
  on.exit(instability:::.cleanup_private_workspace(workspace), add = TRUE)
  files <- instability:::.authoritative_output_files()

  for (filename in unname(files)) {
    data.table::fwrite(
      data.table::data.table(value = filename),
      file.path(workspace, filename),
      sep = "\t"
    )
  }

  result <- instability:::.read_final_outputs(workspace)

  expect_type(result, "list")
  expect_false(inherits(result, "instability"))
  expect_identical(names(result), names(files))
  expect_true(all(vapply(result, data.table::is.data.table, logical(1))))
})

test_that("private workspace is unique and safely cleaned", {
  first <- instability:::.new_private_workspace()
  second <- instability:::.new_private_workspace()
  on.exit(instability:::.cleanup_private_workspace(first), add = TRUE)
  on.exit(instability:::.cleanup_private_workspace(second), add = TRUE)

  expect_false(identical(first, second))
  expect_true(dir.exists(first))
  expect_true(dir.exists(second))

  instability:::.cleanup_private_workspace(first)
  instability:::.cleanup_private_workspace(second)

  expect_false(dir.exists(first))
  expect_false(dir.exists(second))
})
