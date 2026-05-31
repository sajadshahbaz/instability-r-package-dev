#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

infile <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")
outdir <- file.path(nm2_root, "04_results", "figure3_uv_quadrant_locked_v4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clean_numeric <- function(x) {
  as.numeric(gsub("[^0-9eE+.-]", "", as.character(x)))
}

dt <- fread(infile)

if ("feature" %in% names(dt) && !("gene" %in% names(dt))) {
  setnames(dt, "feature", "gene")
}

needed <- c("gene", "delta_pool", "instability_score")
missing_cols <- setdiff(needed, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

for (cc in intersect(c("delta_pool", "instability_score", "hetero_gap", "discordance"), names(dt))) {
  dt[[cc]] <- clean_numeric(dt[[cc]])
}

dt <- copy(dt)
dt[, abs_delta_pool := abs(delta_pool)]

# Locked thresholds
x_cut <- quantile(dt$abs_delta_pool, 0.25, na.rm = TRUE)
y_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)

dt[, quadrant := fifelse(
  abs_delta_pool <= x_cut & instability_score >= y_cut, "masked_by_pooling",
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

# Rank masked genes and keep top 3 for zoom labeling
dt[, masked_rank_score := instability_score - abs_delta_pool]
label_dt <- dt[quadrant == "masked_by_pooling"][
  order(-masked_rank_score, -instability_score, abs_delta_pool)
][1:min(.N, 3)]

x_max <- max(dt$abs_delta_pool, na.rm = TRUE) * 1.03
y_max <- max(dt$instability_score, na.rm = TRUE) * 1.03

base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    plot.margin = margin(8, 8, 8, 8)
  )

# -------- Main panel --------
p_main <- ggplot(dt, aes(x = abs_delta_pool, y = instability_score)) +
  geom_point(aes(color = quadrant), alpha = 0.30, size = 0.9, stroke = 0) +
  geom_vline(xintercept = x_cut, linetype = "dashed", linewidth = 0.5, color = "grey35") +
  geom_hline(xintercept = NULL, yintercept = y_cut, linetype = "dashed", linewidth = 0.5, color = "grey35") +
  annotate(
    "text",
    x = x_cut + (x_max - x_cut) * 0.62,
    y = y_cut + (y_max - y_cut) * 0.78,
    label = "Heterogeneous response",
    color = "grey20",
    size = 4
  ) +
  annotate(
    "text",
    x = x_cut + (x_max - x_cut) * 0.58,
    y = y_cut * 0.34,
    label = "Coherent response",
    color = "grey20",
    size = 4
  ) +
  annotate(
    "text",
    x = x_cut * 1.10,
    y = y_cut * 0.30,
    label = "Background / weak signal",
    color = "grey25",
    size = 3.5,
    hjust = 0
  ) +
  annotate(
    "curve",
    x = x_cut * 1.12,
    y = y_cut + (y_max - y_cut) * 0.93,
    xend = x_cut * 0.50,
    yend = y_cut + (y_max - y_cut) * 0.60,
    curvature = 0.15,
    linewidth = 0.5,
    color = "#B22222",
    arrow = arrow(length = unit(0.10, "inches"), type = "closed")
  ) +
  annotate(
    "text",
    x = x_cut * 1.72,
    y = y_cut + (y_max - y_cut) * 0.97,
    label = "Genes masked by pooling",
    color = "#B22222",
    size = 4.2,
    fontface = "bold",
    hjust = 0
  ) +
    scale_color_manual(values = quad_colors, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max), expand = FALSE) +
  labs(
    title = "UV exposure: pooled effect vs instability",
    subtitle = paste0(
      "Dashed lines mark |pooled effect| Q25 = ", round(x_cut, 3),
      " and instability Q90 = ", round(y_cut, 3)
    ),
    x = "|Pooled log2 fold change|",
    y = "Instability score"
  ) +
  base_theme

# -------- Zoom panel bounds: ONLY masked quadrant --------
x_zoom_min <- 0
x_zoom_max <- x_cut
y_zoom_min <- y_cut
y_zoom_max <- max(dt$instability_score, na.rm = TRUE)

# Keep only labels inside the zoom
label_dt <- label_dt[
  abs_delta_pool >= x_zoom_min & abs_delta_pool <= x_zoom_max &
    instability_score >= y_zoom_min & instability_score <= y_zoom_max
]

# Put labels above the zoom box
if (nrow(label_dt) > 0) {
  setorder(label_dt, abs_delta_pool)
  label_dt[, label_x := seq(0.15, 0.85, length.out = .N) * (x_zoom_max - x_zoom_min) + x_zoom_min]
  label_dt[, label_y := y_zoom_max + (y_zoom_max - y_zoom_min) * 0.08]
} else {
  label_dt[, label_x := numeric()]
  label_dt[, label_y := numeric()]
}

zoom_dt <- dt[
  abs_delta_pool >= x_zoom_min & abs_delta_pool <= x_zoom_max &
    instability_score >= y_zoom_min & instability_score <= y_zoom_max
]

# -------- Zoom panel --------
p_zoom <- ggplot(zoom_dt, aes(x = abs_delta_pool, y = instability_score)) +
  geom_point(
    data = zoom_dt[quadrant != "masked_by_pooling"],
    aes(color = quadrant),
    alpha = 0.30, size = 0.85, stroke = 0
  ) +
  geom_point(
    data = zoom_dt[quadrant == "masked_by_pooling"],
    aes(color = quadrant),
    alpha = 0.55, size = 1.15, stroke = 0
  ) +
  geom_vline(xintercept = x_cut, linetype = "dashed", linewidth = 0.4, color = "grey35") +
  geom_hline(yintercept = y_cut, linetype = "dashed", linewidth = 0.4, color = "grey35") +
    scale_color_manual(values = quad_colors, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(
    xlim = c(x_zoom_min, x_zoom_max),
    ylim = c(y_zoom_min, y_zoom_max + (y_zoom_max - y_zoom_min) * 0.12),
    expand = FALSE
  ) +
  labs(
    title = "Masked region (zoom)",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey94", linewidth = 0.25),
    panel.border = element_rect(linewidth = 0.6, colour = "grey25"),
    plot.margin = margin(34, 8, 8, 8)
  )

# Add labels/arrows only if genes exist
if (nrow(label_dt) > 0) {
  p_zoom <- p_zoom +
    geom_point(
      data = label_dt,
      aes(x = abs_delta_pool, y = instability_score),
      inherit.aes = FALSE,
      color = "#D55E00",
      size = 2.4
    ) +
    geom_text(
      data = label_dt,
      aes(x = label_x, y = label_y, label = gene),
      inherit.aes = FALSE,
      color = "#B22222",
      size = 2.1,
      vjust = 0
    ) +
    geom_segment(
      data = label_dt,
      aes(
        x = label_x,
        y = label_y - 0.01,
        xend = abs_delta_pool,
        yend = instability_score + 0.01
      ),
      inherit.aes = FALSE,
      linewidth = 0.25,
      color = "grey25",
      arrow = arrow(length = unit(0.07, "inches"), type = "closed")
    )
}

# -------- Assemble --------
final_plot <- plot_grid(
  p_main,
  p_zoom,
  nrow = 1,
  rel_widths = c(3.2, 1.08),
  align = "h"
)

png_out <- file.path(outdir, "figure3a_uv_quadrant_locked_v5.png")
pdf_out <- file.path(outdir, "figure3a_uv_quadrant_locked_v5.pdf")

try(
  ggsave(
    png_out,
    final_plot,
    width = 9.0,
    height = 5.8,
    dpi = 400,
    device = "png"
  ),
  silent = TRUE
)

ggsave(
  pdf_out,
  final_plot,
  width = 9.0,
  height = 5.8,
  device = "pdf"
)

if (!file.exists(pdf_out)) {
  stop("PDF figure export failed: ", pdf_out)
}

summary_dt <- data.table(
  condition = "uv",
  x_cut_q25_abs_pooled = round(x_cut, 4),
  y_cut_q90_instability = round(y_cut, 4),
  n_total = nrow(dt),
  n_masked = sum(dt$quadrant == "masked_by_pooling", na.rm = TRUE),
  frac_masked = round(mean(dt$quadrant == "masked_by_pooling", na.rm = TRUE), 4)
)

fwrite(
  summary_dt,
  file.path(outdir, "figure3a_uv_quadrant_locked_v5_summary.tsv"),
  sep = "\t"
)

fwrite(
  label_dt[, .(gene, abs_delta_pool, instability_score, masked_rank_score)],
  file.path(outdir, "figure3a_uv_zoom_labeled_genes.tsv"),
  sep = "\t"
)

cat("Figure 3a UV quadrant locked v5 generated.\n")
print(summary_dt)
cat("\nZoom-labeled genes:\n")
print(label_dt[, .(gene, abs_delta_pool, instability_score, masked_rank_score)])
