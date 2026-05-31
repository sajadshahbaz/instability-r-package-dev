#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ===============================
# CONFIG
# ===============================
nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

in_file <- file.path(
  nm2_root,
  "04_results/figure5b_topN_overlap/figure5b_overlap_robustness.tsv"
)

out_dir <- file.path(nm2_root, "04_results/figure5c_selection_stability_FINAL")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) stop("Missing input file: ", in_file)

# ===============================
# LOAD DATA
# ===============================
dt <- fread(in_file)

if (!all(c("top_n", "jaccard") %in% names(dt))) {
  stop("Required columns missing (top_n, jaccard)")
}

# ===============================
# PLOT (LOCKED DESIGN)
# ===============================
p <- ggplot(dt, aes(x = top_n, y = jaccard)) +
  geom_line(size = 1.2, color = "#1f78b4") +
  geom_point(size = 2, color = "#1f78b4") +
  
  # 🔹 critical reference line (reviewer anchor)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    title = "Figure 5C. Minimal overlap between instability and pooled selection across\n Top-N thresholds",
    subtitle = "Overlap remains near zero across thresholds, confirming persistent divergence",
    x = "Top-N",
    y = "Jaccard index"
  ) +
  
  ylim(0, 1)

# ===============================
# SAVE
# ===============================
try(
  ggsave(
    file.path(out_dir, "Figure5C_selection_stability_FINAL.png"),
    p,
    width = 9,
    height = 5,
    dpi = 300,
    device = "png"
  ),
  silent = TRUE
)

ggsave(
  file.path(out_dir, "Figure5C_selection_stability_FINAL.pdf"),
  p,
  width = 8,
  height = 5,
  device = "pdf"
)

if (!file.exists(file.path(out_dir, "Figure5C_selection_stability_FINAL.pdf"))) {
  stop("PDF export failed for Figure 5C selection-stability panel")
}

# ===============================
# OUTPUT
# ===============================
cat("Figure 5C FINAL generated:\n", out_dir, "\n")
