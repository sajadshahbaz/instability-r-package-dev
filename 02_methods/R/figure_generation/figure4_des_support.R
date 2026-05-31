#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

infile <- file.path(nm2_root, "04_results", "tables", "des_pool_vs_strat_all.tsv")
outdir <- file.path(nm2_root, "04_results", "figure4a_des_support")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clean_numeric <- function(x) {
  as.numeric(gsub("[^0-9eE+.-]", "", as.character(x)))
}

dt <- fread(infile)

if ("feature" %in% names(dt) && !("gene" %in% names(dt))) {
  setnames(dt, "feature", "gene")
}

needed <- c("gene", "delta_pool", "instability_score", "hetero_gap")
missing_cols <- setdiff(needed, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

for (cc in intersect(c("delta_pool", "instability_score", "hetero_gap", "discordance", "var_proxy"), names(dt))) {
  dt[[cc]] <- clean_numeric(dt[[cc]])
}

dt <- dt[
  is.finite(delta_pool) &
  is.finite(instability_score) &
  is.finite(hetero_gap)
]

dt <- copy(dt)
dt[, abs_delta_pool := abs(delta_pool)]

# Locked thresholds
x_cut <- quantile(dt$abs_delta_pool, 0.25, na.rm = TRUE)
y_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)

dt[, masked := abs_delta_pool <= x_cut & instability_score >= y_cut]
dt[, quadrant := fifelse(
  masked, "masked_by_pooling",
  fifelse(
    abs_delta_pool > x_cut & instability_score < y_cut, "coherent_response",
    fifelse(
      abs_delta_pool > x_cut & instability_score >= y_cut, "complex_heterogeneity",
      "background_weak"
    )
  )
)]

quad_colors <- c(
  "background_weak"       = "grey62",
  "coherent_response"     = "grey45",
  "complex_heterogeneity" = "grey32",
  "masked_by_pooling"     = "#D55E00"
)

# tiers for right panel
inst_top10_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)
inst_top5_cut  <- quantile(dt$instability_score, 0.95, na.rm = TRUE)

tier_dt <- rbindlist(list(
  dt[, .(gene, instability_score, hetero_gap, tier = "Background")],
  dt[instability_score >= inst_top10_cut, .(gene, instability_score, hetero_gap, tier = "Top 10%")],
  dt[instability_score >= inst_top5_cut,  .(gene, instability_score, hetero_gap, tier = "Top 5%")]
), use.names = TRUE)

tier_dt[, tier := factor(tier, levels = c("Background", "Top 10%", "Top 5%"))]

tier_summary <- tier_dt[, .(
  n = .N,
  median_hetero_gap = median(hetero_gap, na.rm = TRUE),
  q1 = quantile(hetero_gap, 0.25, na.rm = TRUE),
  q3 = quantile(hetero_gap, 0.75, na.rm = TRUE)
), by = tier]

fwrite(tier_summary, file.path(outdir, "figureS4a_des_tier_summary.tsv"), sep = "\t")

# display limits
x_display_max <- quantile(dt$abs_delta_pool, 0.995, na.rm = TRUE)
y_display_max <- quantile(dt$instability_score, 0.995, na.rm = TRUE)

base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

# Left panel
p_left <- ggplot(dt, aes(x = abs_delta_pool, y = instability_score)) +
  geom_point(aes(color = quadrant), alpha = 0.35, size = 0.9, stroke = 0) +
  scale_color_manual(values = quad_colors, guide = "none") +
  geom_vline(xintercept = x_cut, linetype = "dashed", linewidth = 0.45, color = "grey35") +
  geom_hline(yintercept = y_cut, linetype = "dashed", linewidth = 0.45, color = "grey35") +
  annotate(
    "text",
    x = x_display_max * 0.97,
    y = y_display_max * 0.96,
    label = "DES genes with low pooled effect\ncan still show high instability",
    hjust = 1,
    vjust = 1,
    color = "#B22222",
    size = 3.5,
    fontface = "bold"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(
    xlim = c(0, x_display_max),
    ylim = c(0, y_display_max),
    expand = FALSE
  ) +
  labs(
    title = "DES: pooled effect vs instability",
    x = "|Pooled log2 fold change|",
    y = "Instability score"
  ) +
  base_theme

# Right panel
p_right <- ggplot(tier_dt, aes(x = tier, y = hetero_gap)) +
  geom_boxplot(
    width = 0.62,
    outlier.alpha = 0.10,
    fill = "white",
    color = "grey20",
    linewidth = 0.6
  ) +
  labs(
    title = "High-instability DES genes\n show stronger effect divergence",
    x = NULL,
    y = "Effect divergence"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

final_plot <- plot_grid(
  p_left,
  p_right,
  nrow = 1,
  rel_widths = c(1.65, 1.0),
  align = "h"
)

try(
try(
    ggsave(
      file.path(outdir, "figureS4a_des_support_locked.png"),
      final_plot,
      width = 9.4,
      height = 4.9,
      dpi = 400
    )
,
  silent = TRUE
)
,
  silent = TRUE
)

ggsave(
  file.path(outdir, "figureS4a_des_support_locked.pdf"),
  final_plot,
  width = 9.4,
  height = 4.9,
  device = "pdf"
)

summary_dt <- data.table(
  condition = "des",
  n_total = nrow(dt),
  n_masked = sum(dt$masked, na.rm = TRUE),
  frac_masked = round(mean(dt$masked, na.rm = TRUE), 4),
  x_cut_q25_abs_pooled = round(x_cut, 4),
  y_cut_q90_instability = round(y_cut, 4),
  inst_top10_cut = round(inst_top10_cut, 4),
  inst_top5_cut = round(inst_top5_cut, 4)
)

fwrite(summary_dt, file.path(outdir, "figureS4a_des_summary.tsv"), sep = "\t")

cat("Figure S4A (DES support) generated.\n\n")
print(summary_dt)

cat("\nTier summary:\n")
print(tier_summary)
