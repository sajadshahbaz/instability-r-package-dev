#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

in_file <- file.path(
  nm2_root,
  "04_results", "synthetic_validation_v2", "tables", "synthetic_scenario_summary.tsv"
)

out_dir <- file.path(nm2_root, "04_results", "figure6_variance_disentanglement_FINAL")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) stop("Missing input: ", in_file)

dt <- fread(in_file)
setnames(dt, tolower(names(dt)))

required <- c("scenario", "median_instability", "median_var_proxy", "median_abs_delta_pool")
missing <- setdiff(required, names(dt))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

# keep only the core scenarios for the variance-disentanglement proof
dt_main <- copy(dt)[scenario %in% c("coherent", "opposing", "phase")]

# aggregate across duplicated rows if present (e.g. subgroup/timecourse or noise strata)
dt_main <- dt_main[, .(
  median_instability = median(median_instability, na.rm = TRUE),
  median_var_proxy = median(median_var_proxy, na.rm = TRUE),
  median_abs_delta_pool = median(median_abs_delta_pool, na.rm = TRUE)
), by = scenario]

dt_main[, scenario := factor(scenario, levels = c("coherent", "opposing", "phase"))]

# ----------------------------
# Panel A — instability by scenario
# ----------------------------
p1 <- ggplot(dt_main, aes(x = scenario, y = median_instability, fill = scenario)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Instability rises under opposing internal structure",
    x = "Scenario",
    y = "Median instability"
  )

# ----------------------------
# Panel B — variance vs instability
# ----------------------------
p2 <- ggplot(dt_main, aes(
  x = median_var_proxy,
  y = median_instability,
  color = scenario,
  label = scenario
)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text(nudge_y = 0.15, size = 2.6, show.legend = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Instability is not explained by variance alone",
    x = "Median variance proxy",
    y = "Median instability"
  )

final_plot <- plot_grid(
  p1, p2,
  ncol = 2,
  rel_widths = c(1.05, 1.0)
)

try(
try(
    ggsave(
      file.path(out_dir, "Figure6_variance_disentanglement_FINAL.png"),
      final_plot,
      width = 10,
      height = 4.6,
      dpi = 300
    )
,
  silent = TRUE
)
,
  silent = TRUE
)

ggsave(
  file.path(out_dir, "Figure6_variance_disentanglement_FINAL.pdf"),
  final_plot,
  width = 10,
  height = 4.6
)

fwrite(dt_main, file.path(out_dir, "figure6_source_table.tsv"), sep = "\t")

cat("Figure 6 FINAL generated:\n", out_dir, "\n\n")
print(dt_main)
