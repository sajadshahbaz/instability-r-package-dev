#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

infile <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")
outdir <- file.path(nm2_root, "04_results", "figure5a_rank_displacement")
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

for (cc in intersect(
  c("delta_pool", "instability_score", "hetero_gap", "discordance", "var_proxy"),
  names(dt)
)) {
  dt[[cc]] <- clean_numeric(dt[[cc]])
}

dt <- dt[
  is.finite(delta_pool) &
    is.finite(instability_score)
]

dt <- copy(dt)
dt[, abs_delta_pool := abs(delta_pool)]

# ----------------------------
# Locked thresholds from Fig 3A
# ----------------------------
x_cut <- quantile(dt$abs_delta_pool, 0.25, na.rm = TRUE)
y_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)

dt[, masked := abs_delta_pool <= x_cut & instability_score >= y_cut]

# ----------------------------
# Rank genes
# Rank 1 = strongest
# ----------------------------
dt[, pooled_rank := frank(-abs_delta_pool, ties.method = "average")]
dt[, instability_rank := frank(-instability_score, ties.method = "average")]

# Positive displacement = prioritized more by instability than pooled effect
dt[, rank_displacement := pooled_rank - instability_rank]

# Choose representative masked genes for labels
label_dt <- dt[masked == TRUE][
  order(instability_rank)
][seq(1, .N, length.out = min(.N, 3))]

rank_max <- max(c(dt$pooled_rank, dt$instability_rank), na.rm = TRUE)

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(dt, aes(x = pooled_rank, y = instability_rank)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_point(
    data = dt[masked == FALSE],
    color = "grey75",
    alpha = 0.12,
    size = 0.7,
    stroke = 0
  ) +
  geom_point(
    data = dt[masked == TRUE],
    color = "#D55E00",
    alpha = 0.65,
    size = 1.7,
    stroke = 0
  ) +
  geom_text(
    data = label_dt,
    aes(label = gene),
    color = "black",
    fontface = "bold",
    size = 2.8,
    hjust = 0,
    nudge_x = rank_max * 0.008,
    check_overlap = TRUE
  ) +
  annotate(
    "text",
    x = rank_max * 0.65,
    y = rank_max * 0.75,
    label = "Masked genes are prioritized by instability,\nbut not by pooled effect",
    hjust = 0,
    vjust = 1,
    size = 3.5,
    color = "black",
    fontface = "bold"
  ) +
  scale_x_reverse(expand = c(0.02, 0.02)) +
  scale_y_reverse(expand = c(0.02, 0.02)) +
  labs(
    title = "Instability prioritizes genes missed by pooled analysis",
    subtitle = "Rank 1 indicates strongest signal within each prioritization scheme",
    x = "Rank by pooled |log2 fold change|",
    y = "Rank by instability"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    plot.margin = margin(8, 8, 8, 8)
  )

# ----------------------------
# Save figure
# ----------------------------
png_out <- file.path(outdir, "figure5a_rank_displacement_locked.png")
pdf_out <- file.path(outdir, "figure5a_rank_displacement_locked.pdf")

ggsave(
  filename = png_out,
  plot = p,
  width = 6.8,
  height = 5.6,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = pdf_out,
  plot = p,
  width = 6.8,
  height = 5.6,
  device = "pdf",
  bg = "white"
)

# ----------------------------
# Output tables
# ----------------------------
rank_table <- dt[, .(
  gene,
  abs_delta_pool,
  instability_score,
  pooled_rank,
  instability_rank,
  rank_displacement,
  masked
)][order(instability_rank)]

rank_table_out <- file.path(outdir, "figure5a_rank_table.tsv")
label_out <- file.path(outdir, "figure5a_labeled_genes.tsv")
summary_out <- file.path(outdir, "figure5a_summary.tsv")

fwrite(rank_table, rank_table_out, sep = "\t")

fwrite(
  label_dt[, .(
    gene,
    abs_delta_pool,
    instability_score,
    pooled_rank,
    instability_rank,
    rank_displacement
  )],
  label_out,
  sep = "\t"
)

summary_dt <- data.table(
  n_total = nrow(dt),
  n_masked = sum(dt$masked, na.rm = TRUE),
  x_cut_q25_abs_pooled = round(x_cut, 4),
  y_cut_q90_instability = round(y_cut, 4),
  max_rank = rank_max
)

fwrite(summary_dt, summary_out, sep = "\t")

cat("Figure 5A generated.\n")
cat("PNG:", png_out, "\n")
cat("PDF:", pdf_out, "\n")
cat("Rank table:", rank_table_out, "\n")
cat("Labeled genes:", label_out, "\n")
cat("Summary:", summary_out, "\n\n")

print(summary_dt)

cat("\nTop labeled masked genes:\n")
print(label_dt[, .(
  gene,
  abs_delta_pool,
  instability_score,
  pooled_rank,
  instability_rank,
  rank_displacement
)])
