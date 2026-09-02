#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: compare_original_uv_des.R REPOSITORY_ROOT WORKSPACE", call. = FALSE)
}

repo <- normalizePath(args[[1]], mustWork = TRUE)
workspace <- normalizePath(args[[2]], mustWork = TRUE)
source(file.path(repo, "oracle", "uv_des_expected_schema.R"), local = TRUE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Required package is unavailable: data.table", call. = FALSE)
}

reference_dir <- file.path(repo, "04_results", "tables")
generated_dir <- file.path(workspace, "04_results", "tables")

fail <- function(...) stop(sprintf(...), call. = FALSE)

read_exact <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    fail("Missing or empty table: %s", path)
  }
  data.table::fread(path, na.strings = "NA", check.names = FALSE)
}

numeric_masks <- function(x) list(
  na = is.na(x) & !is.nan(x),
  nan = is.nan(x),
  pos_inf = !is.na(x) & x == Inf,
  neg_inf = !is.na(x) & x == -Inf
)

compare_column <- function(reference, generated, filename, column) {
  if (!identical(typeof(reference), typeof(generated))) {
    fail("%s: type differs for %s (%s vs %s)", filename, column,
         typeof(reference), typeof(generated))
  }
  if (!identical(class(reference), class(generated))) {
    fail("%s: class differs for %s", filename, column)
  }
  if (!identical(length(reference), length(generated))) {
    fail("%s: length differs for %s", filename, column)
  }
  if (is.numeric(reference)) {
    rmasks <- numeric_masks(reference)
    gmasks <- numeric_masks(generated)
    for (mask_name in names(rmasks)) {
      if (!identical(rmasks[[mask_name]], gmasks[[mask_name]])) {
        fail("%s: %s mask differs for %s", filename, mask_name, column)
      }
    }
  }
  if (!identical(reference, generated)) {
    differing <- which(!(reference == generated) | xor(is.na(reference), is.na(generated)))
    first <- if (length(differing)) differing[[1]] else NA_integer_
    fail("%s: exact values differ for %s; first differing row: %s",
         filename, column, first)
  }
}

tables <- list()
for (filename in uv_des_expected$scientific_files) {
  reference <- read_exact(file.path(reference_dir, filename))
  generated <- read_exact(file.path(generated_dir, filename))
  expected_dim <- unname(uv_des_expected$dimensions[[filename]])
  expected_names <- uv_des_schema_for(filename)

  if (!identical(dim(reference), expected_dim)) {
    fail("%s: frozen reference dimensions differ from expected", filename)
  }
  if (!identical(dim(generated), expected_dim)) {
    fail("%s: generated dimensions differ from expected", filename)
  }
  if (!identical(names(reference), expected_names)) {
    fail("%s: frozen reference schema differs from expected", filename)
  }
  if (!identical(names(generated), expected_names)) {
    fail("%s: generated schema differs from expected", filename)
  }
  for (column in expected_names) {
    compare_column(reference[[column]], generated[[column]], filename, column)
  }
  tables[[filename]] <- list(reference = reference, generated = generated)
  message("PARSED EXACT: ", filename)
}

for (condition in c("uv", "des")) {
  all_name <- paste0(condition, "_pool_vs_strat_all.tsv")
  top_name <- paste0(condition, "_pool_vs_strat_top150.tsv")
  for (version in c("reference", "generated")) {
    all_table <- tables[[all_name]][[version]]
    top_table <- tables[[top_name]][[version]]
    calculated <- all_table[order(-instability_score)][seq_len(150L)]
    if (!identical(calculated$feature, top_table$feature)) {
      fail("%s: %s Top-150 membership/order does not follow the original ranking", condition, version)
    }
  }

  for (score in c("hetero_gap", "discordance", "instability_score")) {
    reference <- tables[[all_name]]$reference
    generated <- tables[[all_name]]$generated
    for (k in uv_des_expected$k_values) {
      ref_members <- reference[order(-get(score))][seq_len(min(k, nrow(reference))), feature]
      gen_members <- generated[order(-get(score))][seq_len(min(k, nrow(generated))), feature]
      if (!identical(ref_members, gen_members)) {
        fail("%s: Top-%d membership/order differs for %s", condition, k, score)
      }
    }
  }
}

tail <- tables[["tail_jaccard_metrics.tsv"]]$generated
expected_tail_conditions <- rep(uv_des_expected$condition_order,
                                each = length(uv_des_expected$k_values))
expected_tail_k <- rep(uv_des_expected$k_values,
                       times = length(uv_des_expected$condition_order))
if (!identical(tail$condition, expected_tail_conditions) ||
    !identical(as.integer(tail$K), expected_tail_k)) {
  fail("tail_jaccard_metrics.tsv: condition/K order differs")
}

rank_table <- tables[["rank_displacement_pool_vs_strat.tsv"]]$generated
if (!identical(unique(rank_table$condition), uv_des_expected$condition_order)) {
  fail("rank_displacement_pool_vs_strat.tsv: condition order differs")
}
for (condition_value in uv_des_expected$condition_order) {
  block <- rank_table[rank_table[["condition"]] == condition_value]
  if (!identical(sort(block$rank_strat), seq_len(nrow(block))) ||
      !identical(sort(block$rank_pool), seq_len(nrow(block)))) {
    fail("rank_displacement_pool_vs_strat.tsv: invalid ranks for %s", condition_value)
  }
  if (!identical(block$abs_rank_diff, abs(block$rank_strat - block$rank_pool)) ||
      !identical(block$signed_rank_diff, block$rank_pool - block$rank_strat)) {
    fail("rank_displacement_pool_vs_strat.tsv: inconsistent displacement for %s", condition_value)
  }
}

for (filename in c("pool_vs_strat_metrics.tsv", "tail_headline_summary.tsv",
                   "core_instability_signal_table.tsv")) {
  conditions <- tables[[filename]]$generated$condition
  if (!identical(conditions, uv_des_expected$condition_order)) {
    fail("%s: final condition order differs", filename)
  }
}

message("ALL PARSED SCIENTIFIC COMPARISONS PASSED EXACTLY")
