#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

infile <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")
outdir <- file.path(nm2_root, "04_results", "figure3d_uv_gene_examples_locked")
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

# Detect subgroup columns robustly
cand_short <- c("delta_1", "delta_short", "short_delta", "logFC_short", "delta_A")
cand_long  <- c("delta_2", "delta_long", "long_delta", "logFC_long", "delta_B")

short_col <- cand_short[cand_short %in% names(dt)][1]
long_col  <- cand_long[cand_long %in% names(dt)][1]

if (is.na(short_col) || is.na(long_col)) {
  stop("Subgroup columns not found. Available columns:\n", paste(names(dt), collapse = ", "))
}

for (cc in unique(c("delta_pool", "instability_score", short_col, long_col, "hetero_gap", "discordance"))) {
  if (cc %in% names(dt)) dt[[cc]] <- clean_numeric(dt[[cc]])
}

dt <- dt[
  is.finite(delta_pool) &
  is.finite(instability_score) &
  is.finite(get(short_col)) &
  is.finite(get(long_col))
]

dt[, abs_delta_pool := abs(delta_pool)]

# Locked thresholds from Fig 3A
x_cut <- quantile(dt$abs_delta_pool, 0.25, na.rm = TRUE)
y_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)

dt[, masked := abs_delta_pool <= x_cut & instability_score >= y_cut]

# Strong collapse selection:
# prioritize masked genes with sign flip, strong subgroup opposition, weak pooled effect
dt[, opposition := abs(get(short_col) - get(long_col))]
dt[, sign_flip := sign(get(short_col)) != sign(get(long_col))]
dt[, rank_score := instability_score + opposition - abs_delta_pool]

cand <- dt[masked == TRUE & sign_flip == TRUE][order(-rank_score)]

label_dt <- unique(cand[, .(
  gene,
  short = get(short_col),
  pooled = delta_pool,
  long = get(long_col),
  instability_score,
  opposition,
  abs_delta_pool
)], by = "gene")[1:min(.N, 4)]

if (nrow(label_dt) < 3) stop("Not enough strong collapse genes found.")

fwrite(label_dt, file.path(outdir, "figure3d_selected_genes.tsv"), sep = "\t")

# Reshape for plotting
plot_dt <- melt(
  label_dt,
  id.vars = c("gene", "instability_score", "opposition", "abs_delta_pool"),
  measure.vars = c("short", "pooled", "long"),
  variable.name = "group",
  value.name = "effect"
)

plot_dt[, group := factor(group, levels = c("short", "pooled", "long"))]
plot_dt[, gene := factor(gene, levels = label_dt$gene)]

y_lim <- max(abs(plot_dt$effect), na.rm = TRUE) * 1.20

p <- ggplot(plot_dt, aes(x = group, y = effect, group = gene)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(
    data = plot_dt[group != "pooled"],
    color = "black",
    size = 2.6,
    stroke = 0
  ) +
  geom_point(
    data = plot_dt[group == "pooled"],
    color = "#D55E00",
    size = 4.2,
    stroke = 0
  ) +
  facet_wrap(~ gene, nrow = 1) +
  coord_cartesian(ylim = c(-y_lim, y_lim), expand = FALSE) +
  labs(
    title = "Representative UV genes reveal collapse under pooling",
    subtitle = "Opposing subgroup effects collapse to near-zero under pooling",
    x = NULL,
    y = "Effect size"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

try(
try(
    ggsave(
      file.path(outdir, "figure3d_uv_gene_examples_locked.png"),
      p,
      width = 11.5,
      height = 4.0,
      dpi = 400
    )
,
  silent = TRUE
)
,
  silent = TRUE
)

ggsave(
  file.path(outdir, "figure3d_uv_gene_examples_locked.pdf"),
  p,
  width = 11.5,
  height = 4.0,
  device = "pdf"
)

summary_dt <- data.table(
  n_total = nrow(dt),
  n_masked = sum(dt$masked, na.rm = TRUE),
  x_cut_q25_abs_pooled = round(x_cut, 4),
  y_cut_q90_instability = round(y_cut, 4),
  short_column = short_col,
  long_column = long_col
)

fwrite(summary_dt, file.path(outdir, "figure3d_summary.tsv"), sep = "\t")

cat("Figure 3D LOCKED generated.\n\n")
print(summary_dt)

cat("\nSelected genes:\n")
print(label_dt)
