#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

in_dir  <- file.path(nm2_root, "04_results", "tables")
out_dir <- file.path(nm2_root, "04_results", "gate3_instability_vs_variance")

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

clean_numeric <- function(x){
  as.numeric(gsub("[^0-9eE+.-]", "", as.character(x)))
}

run_condition <- function(cond){

  cat("Processing:", cond, "\n")

  path <- file.path(in_dir, paste0(cond,"_pool_vs_strat_all.tsv"))
  dt <- fread(path)

  numeric_cols <- intersect(c("instability_score","var_proxy","hetero_gap","discordance"), names(dt))
  for(cc in numeric_cols){
    dt[[cc]] <- clean_numeric(dt[[cc]])
  }

  if ("sign_flip" %in% names(dt)) {
    dt[, sign_flip_logical := as.character(sign_flip) %in% c("TRUE","True","true","1")]
  } else if ("discordance" %in% names(dt)) {
    dt[, sign_flip_logical := discordance > 0]
  } else {
    dt[, sign_flip_logical := FALSE]
  }

  spearman_cor <- cor(
    dt$instability_score,
    dt$var_proxy,
    method = "spearman",
    use = "complete.obs"
  )

  p <- ggplot(dt, aes(x = var_proxy, y = instability_score)) +
    geom_point(alpha = 0.25, size = 0.8) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(toupper(cond), ": instability vs variance"),
      subtitle = paste0("Spearman rho = ", round(spearman_cor, 3)),
      x = "variance proxy",
      y = "instability score"
    )

  ggsave(
    file.path(fig_dir, paste0(cond,"_instability_vs_variance.png")),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )

  inst_cut <- quantile(dt$instability_score, 0.90, na.rm=TRUE)
  var_cut  <- quantile(dt$var_proxy, 0.90, na.rm=TRUE)

  dt[, top_instability := instability_score >= inst_cut]
  dt[, top_variance    := var_proxy >= var_cut]

  id_col <- if ("feature" %in% names(dt)) "feature" else if ("gene" %in% names(dt)) "gene" else stop("No feature/gene column found")

  inst_set <- dt[top_instability == TRUE][[id_col]]
  var_set  <- dt[top_variance == TRUE][[id_col]]

  overlap <- length(intersect(inst_set, var_set))
  overlap_fraction <- overlap / length(inst_set)

  inst_flip_rate   <- mean(dt[top_instability == TRUE]$sign_flip_logical, na.rm=TRUE)
  global_flip_rate <- mean(dt$sign_flip_logical, na.rm=TRUE)

  summary <- data.table(
    condition = cond,
    n_features = nrow(dt),
    spearman_instability_vs_variance = round(spearman_cor, 4),
    top10_instability = length(inst_set),
    top10_variance = length(var_set),
    overlap_top_sets = overlap,
    overlap_fraction = round(overlap_fraction, 4),
    global_signflip_rate = round(global_flip_rate, 4),
    instability_signflip_rate = round(inst_flip_rate, 4)
  )

  fwrite(
    summary,
    file.path(tab_dir, paste0(cond,"_instability_vs_variance_summary.tsv")),
    sep="\t"
  )

  return(summary)
}

results <- rbindlist(
  lapply(c("uv","des","gam"), run_condition)
)

fwrite(
  results,
  file.path(tab_dir,"gate3_instability_vs_variance_summary.tsv"),
  sep="\t"
)

cat("\nGate 3 completed.\n")
