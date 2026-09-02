.new_private_workspace <- function() {
  workspace <- tempfile("instability-")
  if (!dir.create(workspace, recursive = FALSE, showWarnings = FALSE)) {
    stop("Unable to create private instability workspace", call. = FALSE)
  }
  normalizePath(workspace, mustWork = TRUE)
}

.authoritative_output_files <- function() {
  c(
    uv_pool_vs_strat_all = "uv_pool_vs_strat_all.tsv",
    uv_pool_vs_strat_top150 = "uv_pool_vs_strat_top150.tsv",
    uv_pool_vs_strat_summary = "uv_pool_vs_strat_summary.tsv",
    des_pool_vs_strat_all = "des_pool_vs_strat_all.tsv",
    des_pool_vs_strat_top150 = "des_pool_vs_strat_top150.tsv",
    des_pool_vs_strat_summary = "des_pool_vs_strat_summary.tsv",
    pool_vs_strat_metrics = "pool_vs_strat_metrics.tsv",
    tail_jaccard_metrics = "tail_jaccard_metrics.tsv",
    rank_displacement_pool_vs_strat = "rank_displacement_pool_vs_strat.tsv",
    tail_headline_summary = "tail_headline_summary.tsv",
    core_instability_signal_table = "core_instability_signal_table.tsv"
  )
}

.write_condition_outputs <- function(
  workspace,
  condition,
  condition_table,
  top150,
  summary
) {
  data.table::fwrite(
    condition_table,
    file.path(workspace, paste0(condition, "_pool_vs_strat_all.tsv")),
    sep = "\t"
  )
  data.table::fwrite(
    top150,
    file.path(workspace, paste0(condition, "_pool_vs_strat_top150.tsv")),
    sep = "\t"
  )
  data.table::fwrite(
    summary,
    file.path(workspace, paste0(condition, "_pool_vs_strat_summary.tsv")),
    sep = "\t"
  )

  invisible(NULL)
}

.read_final_outputs <- function(workspace) {
  files <- .authoritative_output_files()
  result <- lapply(
    unname(files),
    function(filename) data.table::fread(file.path(workspace, filename))
  )
  names(result) <- names(files)
  result
}

.cleanup_private_workspace <- function(workspace) {
  if (!is.character(workspace) || length(workspace) != 1L || !nzchar(workspace)) {
    stop("Invalid private instability workspace", call. = FALSE)
  }
  if (!dir.exists(workspace)) {
    return(invisible(NULL))
  }

  resolved <- normalizePath(workspace, mustWork = TRUE)
  temp_root <- normalizePath(tempdir(), mustWork = TRUE)
  if (!identical(dirname(resolved), temp_root) ||
      !startsWith(basename(resolved), "instability-")) {
    stop("Refusing to remove an invalid private instability workspace", call. = FALSE)
  }

  unlink(resolved, recursive = TRUE, force = FALSE)
  if (dir.exists(resolved)) {
    stop("Unable to remove private instability workspace", call. = FALSE)
  }

  invisible(NULL)
}
