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
outdir <- file.path(nm2_root, "04_results", "figure3c_uv_instability_vs_disagreement")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clean_numeric <- function(x) {
  as.numeric(gsub("[^0-9eE+.-]", "", as.character(x)))
}

# ----------------------------
# Load and clean
# ----------------------------
dt <- fread(infile)

if ("feature" %in% names(dt) && !("gene" %in% names(dt))) {
  setnames(dt, "feature", "gene")
}

needed <- c("gene", "delta_pool", "instability_score", "discordance")
missing_cols <- setdiff(needed, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

for (cc in intersect(c("delta_pool", "instability_score", "discordance"), names(dt))) {
  dt[[cc]] <- clean_numeric(dt[[cc]])
}

dt <- dt[
  is.finite(delta_pool) &
  is.finite(instability_score) &
  is.finite(discordance)
]

dt[, abs_delta_pool := abs(delta_pool)]

# ----------------------------
# Locked thresholds (from Fig 3A)
# ----------------------------
x_cut <- quantile(dt$abs_delta_pool, 0.25, na.rm = TRUE)
y_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)

dt[, masked := abs_delta_pool <= x_cut & instability_score >= y_cut]

# ----------------------------
# Instability tiers
# ----------------------------
inst_top10_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)
inst_top5_cut  <- quantile(dt$instability_score, 0.95, na.rm = TRUE)

tier_dt <- rbindlist(list(
  dt[, .(gene, instability_score, discordance, tier = "background_all")],
  dt[instability_score >= inst_top10_cut, .(gene, instability_score, discordance, tier = "top10_instability")],
  dt[instability_score >= inst_top5_cut,  .(gene, instability_score, discordance, tier = "top5_instability")]
))

tier_dt[, tier := factor(
  tier,
  levels = c("background_all", "top10_instability", "top5_instability"),
  labels = c("Background", "Top 10%", "Top 5%")
)]

# ----------------------------
# Display limits (CRITICAL FIX)
# ----------------------------
x_display_max <- quantile(dt$instability_score, 0.995, na.rm = TRUE)
y_display_max <- quantile(dt$discordance, 0.99, na.rm = TRUE)

# ----------------------------
# Theme
# ----------------------------
base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

# ----------------------------
# LEFT PANEL (FIXED)
# ----------------------------
p_left <- ggplot(dt, aes(x = instability_score, y = discordance)) +

  geom_point(
    data = dt[masked == FALSE],
    color = "grey70",
    alpha = 0.28,
    size = 1.0
  ) +

  geom_point(
    data = dt[masked == TRUE],
    color = "#D55E00",
    alpha = 0.75,
    size = 1.3
  ) +

  annotate(
    "text",
    x = x_display_max * 0.05,
    y = y_display_max * 0.95,
    label = "Masked genes cluster in\nhigh-instability, high-conflict regions",
    hjust = 0,
    vjust = 1,
    color = "#B22222",
    size = 3.5,
    fontface = "bold"
  ) +

  coord_cartesian(
    xlim = c(0, x_display_max),
    ylim = c(0, y_display_max),
    expand = FALSE
  ) +

  labs(
    title = "Instability captures subgroup conflict",
    x = "Instability score",
    y = "Directional disagreement (discordance)"
  ) +

  base_theme

# ----------------------------
# RIGHT PANEL (unchanged, already good)
# ----------------------------
p_right <- ggplot(tier_dt, aes(x = tier, y = discordance)) +

  geom_boxplot(
    width = 0.6,
    fill = "white",
    color = "grey20",
    outlier.alpha = 0.1
  ) +

  labs(
    title = "High instability tiers show\n stronger directional conflict",
    x = NULL,
    y = "discordance"
  ) +

  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

# ----------------------------
# Assemble
# ----------------------------
final_plot <- plot_grid(
  p_left,
  p_right,
  nrow = 1,
  rel_widths = c(1.65, 1.0)
)

# ----------------------------
# Save
# ----------------------------
try(
try(
    ggsave(
      file.path(outdir, "figure3c_uv_instability_vs_disagreement_locked_v2.png"),
      final_plot,
      width = 9.5,
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
  file.path(outdir, "figure3c_uv_instability_vs_disagreement_locked_v2.pdf"),
  final_plot,
  width = 9.5,
  height = 4.9,
  device = "pdf"
)

cat("Figure 3C (FINAL v2) generated.\n")
