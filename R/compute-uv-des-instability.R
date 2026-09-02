compute_uv_des_instability <- function(vst, metadata) {
  prepared <- .prepare_uv_des_inputs(vst, metadata)

  workspace <- .new_private_workspace()
  on.exit(.cleanup_private_workspace(workspace), add = TRUE)

  uv <- .compute_uv_des_condition(prepared$X, prepared$meta, "uv")
  uv_top150 <- .make_top150(uv)
  uv_summary <- .summarize_uv_des_condition(uv)
  .write_condition_outputs(
    workspace,
    "uv",
    uv,
    uv_top150,
    uv_summary
  )

  des <- .compute_uv_des_condition(prepared$X, prepared$meta, "des")
  des_top150 <- .make_top150(des)
  des_summary <- .summarize_uv_des_condition(des)
  .write_condition_outputs(
    workspace,
    "des",
    des,
    des_top150,
    des_summary
  )

  .make_pool_strat_metrics(workspace)
  .make_tail_outputs(workspace)
  .make_core_table(workspace)

  .read_final_outputs(workspace)
}
