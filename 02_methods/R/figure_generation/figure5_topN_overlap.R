#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

infile <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")
outdir <- file.path(nm2_root, "04_results", "figure5b_topN_overlap")
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

needed <- c("gene", "delta_pool", "instability_score")
missing_cols <- setdiff(needed, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

for (cc in intersect(c("delta_pool", "instability_score"), names(dt))) {
  dt[[cc]] <- clean_numeric(dt[[cc]])
}

dt <- dt[
  is.finite(delta_pool) &
  is.finite(instability_score)
]

dt <- unique(dt, by = "gene")
dt[, abs_delta_pool := abs(delta_pool)]

# ----------------------------
# Rank genes
# rank 1 = strongest
# ----------------------------
dt[, pooled_rank := frank(-abs_delta_pool, ties.method = "average")]
dt[, instability_rank := frank(-instability_score, ties.method = "average")]

# ----------------------------
# Main Figure 5B: Top N overlap
# ----------------------------
TOP_N <- 100

pooled_top <- dt[order(pooled_rank)][1:TOP_N, gene]
inst_top   <- dt[order(instability_rank)][1:TOP_N, gene]

overlap <- intersect(pooled_top, inst_top)
pooled_only <- setdiff(pooled_top, inst_top)
inst_only   <- setdiff(inst_top, pooled_top)

summary_dt <- data.table(
  category = c("Overlap", "Instability-only", "Pooled-only"),
  count = c(length(overlap), length(inst_only), length(pooled_only))
)

summary_dt[, category := factor(
  category,
  levels = c("Overlap", "Instability-only", "Pooled-only")
)]

# main plot
p_main <- ggplot(summary_dt, aes(x = category, y = count, fill = category)) +
  geom_col(width = 0.62) +
  geom_text(
    aes(label = count),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "Overlap" = "grey50",
    "Instability-only" = "#D55E00",
    "Pooled-only" = "#0072B2"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(
    title = paste0("Instability and pooled analysis identify\n non-overlapping gene sets"),
    subtitle = "Instability identifies a largely distinct set of genes",
    x = NULL,
    y = "Number of genes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    plot.margin = margin(8, 8, 8, 8)
  )

try(
try(
    ggsave(
      file.path(outdir, "figure5b_topN_overlap_locked.png"),
      p_main,
      width = 5.8,
      height = 5.0,
      dpi = 400
    )
,
  silent = TRUE
)
,
  silent = TRUE
)

ggsave(
  file.path(outdir, "figure5b_topN_overlap_locked.pdf"),
  p_main,
  width = 5.8,
  height = 5.0,
  device = "pdf"
)

# ----------------------------
# Save main tables
# ----------------------------
fwrite(summary_dt, file.path(outdir, "figure5b_summary.tsv"), sep = "\t")

fwrite(data.table(gene = overlap),
       file.path(outdir, "overlap_genes.tsv"), sep = "\t")

fwrite(data.table(gene = inst_only),
       file.path(outdir, "instability_only_genes.tsv"), sep = "\t")

fwrite(data.table(gene = pooled_only),
       file.path(outdir, "pooled_only_genes.tsv"), sep = "\t")

# ----------------------------
# Robustness across Top N
# ----------------------------
top_ns <- c(25, 50, 100, 200, 500)

robust_dt <- rbindlist(lapply(top_ns, function(n) {
  pooled_n <- dt[order(pooled_rank)][1:min(n, .N), gene]
  inst_n   <- dt[order(instability_rank)][1:min(n, .N), gene]

  overlap_n <- length(intersect(pooled_n, inst_n))
  union_n   <- length(union(pooled_n, inst_n))

  data.table(
    top_n = n,
    overlap = overlap_n,
    pooled_only = length(setdiff(pooled_n, inst_n)),
    instability_only = length(setdiff(inst_n, pooled_n)),
    jaccard = ifelse(union_n > 0, overlap_n / union_n, NA_real_)
  )
}))

fwrite(robust_dt, file.path(outdir, "figure5b_overlap_robustness.tsv"), sep = "\t")

# long table for supplementary plot
robust_dt[, `:=`(
  overlap = as.numeric(overlap),
  jaccard = as.numeric(jaccard)
)]

robust_long <- melt(
  robust_dt,
  id.vars = "top_n",
  measure.vars = c("overlap", "jaccard"),
  variable.name = "metric",
  value.name = "value"
)

robust_long[, metric := factor(
  metric,
  levels = c("overlap", "jaccard"),
  labels = c("Overlap count", "Jaccard index")
)]

# supplementary robustness plot
p_sup <- ggplot(robust_long, aes(x = top_n, y = value)) +
  geom_line(linewidth = 0.7, color = "#333333") +
  geom_point(size = 2.4, color = "#D55E00") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = top_ns) +
  labs(
    title = "Set divergence remains stable across Top-N thresholds",
    subtitle = "Overlap between pooled-effect and instability-selected genes remains limited across\n ranking cutoffs",
    x = "Top N cutoff",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    plot.margin = margin(8, 8, 8, 8)
  )

try(
try(
    ggsave(
      file.path(outdir, "figure5b_overlap_robustness_supplementary.png"),
      p_sup,
      width = 6.4,
      height = 6.0,
      dpi = 400
    )
,
  silent = TRUE
)
,
  silent = TRUE
)

ggsave(
  file.path(outdir, "figure5b_overlap_robustness_supplementary.pdf"),
  p_sup,
  width = 6.4,
  height = 6.0,
  device = "pdf"
)

# ----------------------------
# Console summary
# ----------------------------
cat("Figure 5B generated.\n\n")
cat("Main Top-N summary (N =", TOP_N, "):\n")
print(summary_dt)

cat("\nRobustness across Top-N:\n")
print(robust_dt)
