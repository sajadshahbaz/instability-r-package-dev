#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

outdir <- file.path(nm2_root, "04_results", "figure4_mixed_real_and_control_conditions_LOCKED")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clean_numeric <- function(x) {
  as.numeric(gsub("[^0-9eE+.-]", "", as.character(x)))
}

read_table_generic <- function(cond_name, filename, source_type = c("real", "null")) {
  source_type <- match.arg(source_type)
  infile <- file.path(nm2_root, "04_results", "tables", filename)
  if (!file.exists(infile)) stop("Missing file: ", infile)

  dt <- fread(infile)

  if ("feature" %in% names(dt) && !("gene" %in% names(dt))) {
    setnames(dt, "feature", "gene")
  }

  needed <- c("gene", "instability_score", "hetero_gap", "discordance")
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", infile, ": ", paste(missing_cols, collapse = ", "))
  }

  for (cc in intersect(c("instability_score", "hetero_gap", "discordance", "delta_pool", "var_proxy"), names(dt))) {
    dt[[cc]] <- clean_numeric(dt[[cc]])
  }

  dt <- dt[
    is.finite(instability_score) &
      is.finite(hetero_gap) &
      is.finite(discordance)
  ]

  dt[, condition := cond_name]
  dt[, source_type := source_type]
  dt
}

tables_available <- function(fname) {
  file.exists(file.path(nm2_root, "04_results", "tables", fname))
}

objs <- list()

# Structured conditions
objs[["UV"]]  <- read_table_generic("UV",  "uv_pool_vs_strat_all.tsv", "real")
objs[["DES"]] <- read_table_generic("DES", "des_pool_vs_strat_all.tsv", "real")
objs[["GAM"]] <- read_table_generic("GAM", "gam_pool_vs_strat_all.tsv", "real")

# Homogeneous null control
if (tables_available("lt_null_control_all.tsv")) {
  objs[["LT"]] <- read_table_generic("LT", "lt_null_control_all.tsv", "null")
}

all_dt <- rbindlist(objs, use.names = TRUE, fill = TRUE)
present_levels <- names(objs)
all_dt[, condition := factor(condition, levels = present_levels)]

global_top10_cut  <- quantile(all_dt$instability_score, 0.90, na.rm = TRUE)
global_top5_cut   <- quantile(all_dt$instability_score, 0.95, na.rm = TRUE)
global_pooled_q25 <- if ("delta_pool" %in% names(all_dt)) {
  quantile(abs(all_dt$delta_pool), 0.25, na.rm = TRUE)
} else {
  NA_real_
}

tail_summary <- all_dt[, .(
  n_total = .N,
  frac_top10 = mean(instability_score >= global_top10_cut, na.rm = TRUE),
  frac_top5  = mean(instability_score >= global_top5_cut,  na.rm = TRUE)
), by = .(condition, source_type)]

tail_long <- melt(
  tail_summary,
  id.vars = c("condition", "source_type", "n_total"),
  measure.vars = c("frac_top10", "frac_top5"),
  variable.name = "tail",
  value.name = "fraction"
)

tail_long[, tail := factor(
  tail,
  levels = c("frac_top10", "frac_top5"),
  labels = c("Top 10%", "Top 5%")
)]

metric_summary <- all_dt[, .(
  median_instability = median(instability_score, na.rm = TRUE),
  median_hetero_gap  = median(hetero_gap, na.rm = TRUE),
  median_discordance = median(discordance, na.rm = TRUE)
), by = .(condition, source_type)]

metric_long <- melt(
  metric_summary,
  id.vars = c("condition", "source_type"),
  measure.vars = c("median_instability", "median_hetero_gap", "median_discordance"),
  variable.name = "metric",
  value.name = "value"
)

metric_long[, metric := factor(
  metric,
  levels = c("median_instability", "median_hetero_gap", "median_discordance"),
  labels = c("Median instability", "Median hetero_gap", "Median discordance")
)]

masked_summary <- if ("delta_pool" %in% names(all_dt)) {
  all_dt[, .(
    n_total = .N,
    masked_frac = mean(
      abs(delta_pool) < global_pooled_q25 &
        instability_score >= global_top10_cut,
      na.rm = TRUE
    )
  ), by = .(condition, source_type)]
} else {
  data.table(
    condition = factor(present_levels, levels = present_levels),
    source_type = c(rep("real", length(present_levels))),
    n_total = NA_integer_,
    masked_frac = NA_real_
  )
}

fwrite(tail_summary,   file.path(outdir, "figure4_tail_summary.tsv"), sep = "\t")
fwrite(metric_summary, file.path(outdir, "figure4_metric_summary.tsv"), sep = "\t")
fwrite(masked_summary, file.path(outdir, "figure4_masked_summary.tsv"), sep = "\t")

real_fill <- "#D9D9D9"
null_fill <- "#F4F4F4"
accent_tail <- "#D55E00"
accent_mask <- "#B22222"

base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 6)),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "top",
    legend.title = element_blank()
  )

# a. Instability distribution
p4a <- ggplot(all_dt, aes(x = condition, y = instability_score)) +
  geom_boxplot(
    aes(fill = source_type),
    width = 0.62, outlier.alpha = 0.08, color = "grey20", linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(real = real_fill, null = null_fill),
    labels = c(real = "Structured heterogeneity", null = "Homogeneous null")
  ) +
  labs(title = "a. Instability distribution", x = NULL, y = "Instability score") +
  annotate(
    "text",
    x = 4,
    y = max(all_dt$instability_score, na.rm = TRUE) * 0.88,
    label = "Homogeneous data \u2192 low instability",
    size = 3.5, hjust = 1, vjust= -1.5, color = "grey20"
  ) +
  base_theme

# b. Instability tail occupancy
p4b <- ggplot(tail_long, aes(x = condition, y = fraction, group = tail)) +
  geom_line(color = "grey35", linewidth = 0.2) +
  geom_point(aes(shape = tail), size = 2.8, color = accent_tail) +
  scale_shape_manual(values = c(16, 17)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(title = "b. Instability tail occupancy", x = NULL, y = "Fraction of genes") +
  annotate(
    "text",
    x = 3,
    y = max(tail_long$fraction, na.rm = TRUE) * 1.03,
    label = "Structured conditions enrich instability tail",
    size = 3.3, hjust = 0.5, vjust = 0.1, color = "grey20"
  ) +
  base_theme

# c. Metric decomposition
p4c <- ggplot(metric_long, aes(x = condition, y = value, group = metric)) +
  geom_line(color = "grey35", linewidth = 0.2) +
  geom_point(aes(shape = metric), size = 2.8, color = "black") +
  scale_shape_manual(values = c(16, 17, 15), guide = guide_legend(nrow = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(title = "c. Metric decomposition", x = NULL, y = "Median metric value") +
  annotate(
    "text",
    x = 2.2,
    y = max(metric_long$value, na.rm = TRUE) * 0.95,
    label = "Instability arises from distinct structural drivers",
    size = 3.3, hjust = 0.1, vjust = -2, color = "grey20"
  ) +
  base_theme

# d. Masked fraction
p4d <- ggplot(masked_summary, aes(x = condition, y = masked_frac, group = 1)) +
  geom_line(color = "grey35", linewidth = 0.2) +
  geom_point(size = 2.8, color = accent_mask) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(title = "d. Masked-gene fraction", x = NULL, y = "Fraction of genes") +
  annotate(
    "text",
    x = 3,
    y = max(masked_summary$masked_frac, na.rm = TRUE) * 1.03,
    label = "Pooling masks signal only in structured data",
    size = 3.3, hjust = 0.6, vjust = -0.1, color = "grey20"
  ) +
  base_theme

title_grob <- ggdraw() +
  draw_label(
    "Instability separates structured heterogeneity from homogeneous data",
    fontface = "bold",
    size = 15,
    hjust = 0.5,
    x = 0.5
  )

subtitle_grob <- ggdraw() +
  draw_label(
    "UV: conflict | DES: magnitude divergence | GAM: temporal | LT: null",
    size = 11,
    hjust = 0.5,
    x = 0.5
  )

top_row <- plot_grid(p4a, p4b, nrow = 1, rel_widths = c(1.18, 0.98), align = "h")
bottom_row <- plot_grid(p4c, p4d, nrow = 1, rel_widths = c(1.05, 0.95), align = "h")
panels <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1.0, 1.0))

final_plot <- plot_grid(
  title_grob,
  subtitle_grob,
  panels,
  ncol = 1,
  rel_heights = c(0.07, 0.05, 1.0)
)

try(
try(
    ggsave(
      file.path(outdir, "figure4_mixed_real_and_control_conditions_LOCKED.png"),
      final_plot, width = 14.0, height = 9.4, dpi = 400
    )
,
  silent = TRUE
)
,
  silent = TRUE
)

ggsave(
  file.path(outdir, "figure4_mixed_real_and_control_conditions_LOCKED.pdf"),
  final_plot, width = 14.0, height = 9.4, device = "pdf"
)

cat("LOCKED Figure 4D generated.\n\n")
cat("Conditions included:\n")
print(present_levels)

cat("\nGlobal thresholds:\n")
cat("  top10 instability cut =", round(global_top10_cut, 4), "\n")
cat("  top5  instability cut =", round(global_top5_cut, 4), "\n")
if (is.finite(global_pooled_q25)) {
  cat("  pooled |delta_pool| Q25 =", round(global_pooled_q25, 4), "\n")
}

cat("\nTail summary:\n")
print(tail_summary)

cat("\nMetric summary:\n")
print(metric_summary)

cat("\nMasked summary:\n")
print(masked_summary)
