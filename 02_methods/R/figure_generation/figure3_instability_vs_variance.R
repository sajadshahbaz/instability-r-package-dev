#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(broom)
  library(grid)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

infile <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")
outdir <- file.path(nm2_root, "04_results", "figure3b_instability_vs_variance")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clean_numeric <- function(x) {
  as.numeric(gsub("[^0-9eE+.-]", "", as.character(x)))
}

# ----------------------------
# Load and clean data
# ----------------------------
dt <- fread(infile)

if ("feature" %in% names(dt) && !("gene" %in% names(dt))) {
  setnames(dt, "feature", "gene")
}

needed <- c("gene", "delta_pool", "instability_score", "var_proxy")
missing_cols <- setdiff(needed, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

for (cc in intersect(c("delta_pool", "instability_score", "var_proxy", "hetero_gap", "discordance"), names(dt))) {
  dt[[cc]] <- clean_numeric(dt[[cc]])
}

dt <- dt[is.finite(delta_pool) & is.finite(instability_score) & is.finite(var_proxy)]
dt <- copy(dt)
dt[, abs_delta_pool := abs(delta_pool)]

# ----------------------------
# Locked thresholds from Fig 3A
# ----------------------------
x_cut <- quantile(dt$abs_delta_pool, 0.25, na.rm = TRUE)
y_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)

dt[, masked := abs_delta_pool <= x_cut & instability_score >= y_cut]
dt[, masked_label := fifelse(masked, "TRUE", "FALSE")]

# ----------------------------
# Logistic model
# masked ~ instability + variance
# ----------------------------
fit <- glm(masked ~ instability_score + var_proxy, data = dt, family = binomial())

capture.output(summary(fit), file = file.path(outdir, "figure3b_logistic_summary.txt"))

coef_dt <- as.data.table(tidy(fit, conf.int = TRUE))
coef_dt <- coef_dt[term != "(Intercept)"]
coef_dt[, term_clean := fifelse(
  term == "instability_score", "Instability",
  fifelse(term == "var_proxy", "Variance", term)
)]

fwrite(coef_dt, file.path(outdir, "figure3b_logistic_coefficients.tsv"), sep = "\t")

or_dt <- copy(coef_dt)
or_dt[, `:=`(
  odds_ratio = exp(estimate),
  or_low = exp(conf.low),
  or_high = exp(conf.high)
)]
fwrite(or_dt, file.path(outdir, "figure3b_logistic_odds_ratios.tsv"), sep = "\t")

# ----------------------------
# Data-driven display limits
# Keep model on full data, trim only the visible tail
# ----------------------------
x_display_max <- quantile(dt$var_proxy, 0.995, na.rm = TRUE)
x_display_max <- max(x_display_max, 1.5)

y_display_max <- max(dt$instability_score, na.rm = TRUE) * 1.02

# ----------------------------
# Themes
# ----------------------------
base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    plot.margin = margin(8, 8, 8, 8)
  )

# ----------------------------
# LEFT PANEL: scatter
# ----------------------------
p_scatter <- ggplot(dt, aes(x = var_proxy, y = instability_score)) +
  geom_point(
    data = dt[masked == FALSE],
    color = "grey70",
    alpha = 0.35,
    size = 1.0,
    stroke = 0
  ) +
  geom_point(
    data = dt[masked == TRUE],
    color = "#D55E00",
    alpha = 0.60,
    size = 1.25,
    stroke = 0
  ) +
  annotate(
  "text",
  x = x_display_max * 0.75,
  y = y_display_max * 0.92,
  label = "Masked genes concentrate\nat high instability,\nnot high variance",
  hjust = 1,
  vjust = 1,
  size = 3.6,
  color = "#B22222",
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
    title = "Instability captures masking beyond variance",
    x = "Variance proxy",
    y = "Instability score"
  ) +
  base_theme

# ----------------------------
# RIGHT PANEL: coefficient plot
# ----------------------------
p_coef <- ggplot(coef_dt, aes(x = term_clean, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, color = "grey35") +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.10,
    linewidth = 0.6,
    color = "grey20"
  ) +
  geom_point(size = 3.0, color = "#D55E00") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.10))) +
  labs(
    title = "Instability predicts\n masking beyond variance",
    x = NULL,
    y = "Log-odds coefficient"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    plot.margin = margin(8, 8, 8, 8)
  )

# ----------------------------
# Assemble
# ----------------------------
final_plot <- plot_grid(
  p_scatter,
  p_coef,
  nrow = 1,
  rel_widths = c(1.95, 1.0),
  align = "h"
)

try(
try(
    ggsave(
      file.path(outdir, "figure3b_instability_vs_variance_locked.png"),
      final_plot,
      width = 9.2,
      height = 4.8,
      dpi = 400
    )
,
  silent = TRUE
)
,
  silent = TRUE
)

ggsave(
  file.path(outdir, "figure3b_instability_vs_variance_locked.pdf"),
  final_plot,
  width = 9.2,
  height = 4.8,
  device = "pdf"
)

# ----------------------------
# Console summary
# ----------------------------
summary_dt <- data.table(
  condition = "uv",
  n_total = nrow(dt),
  n_masked = sum(dt$masked, na.rm = TRUE),
  frac_masked = round(mean(dt$masked, na.rm = TRUE), 4),
  x_cut_q25_abs_pooled = round(x_cut, 4),
  y_cut_q90_instability = round(y_cut, 4),
  x_display_max = round(x_display_max, 4)
)

cat("Figure 3B generated.\n\n")
print(summary_dt)

cat("\nLogistic coefficients:\n")
print(coef_dt[, .(term_clean, estimate, conf.low, conf.high, p.value)])
