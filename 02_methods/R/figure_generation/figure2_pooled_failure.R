#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(limma)
})

# =========================================================
# Setup
# =========================================================
nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

out_dir <- file.path(nm2_root, "04_results", "figure2b")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

# =========================================================
# Locked design
# =========================================================
n_null     <- 700
n_coherent <- 100
n_opposing <- 100
n_phase    <- 100
n_genes    <- n_null + n_coherent + n_opposing + n_phase

genes <- c(
  paste0("null_", seq_len(n_null)),
  paste0("coherent_", seq_len(n_coherent)),
  paste0("opposing_", seq_len(n_opposing)),
  paste0("phase_", seq_len(n_phase))
)

truth <- data.table(
  gene = genes,
  scenario = c(
    rep("null", n_null),
    rep("coherent", n_coherent),
    rep("opposing", n_opposing),
    rep("phase", n_phase)
  )
)

n_rep <- 3
n_sim <- 100
noise_sd <- 0.50
timepoints <- c("t1", "t2", "t3", "t4", "t5")

# =========================================================
# Limma helpers
# =========================================================
run_pooled_limma_subgroup <- function(expr_mat, sample_info) {
  design <- model.matrix(~ 0 + group, data = sample_info)
  colnames(design) <- sub("^group", "", colnames(design))
  fit <- lmFit(expr_mat, design)
  cont <- makeContrasts(Treatment - Control, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), pooled_logFC = tt$logFC)
}

run_stratified_limma_subgroup <- function(expr_mat, sample_info) {
  design <- model.matrix(~ 0 + label, data = sample_info)
  colnames(design) <- sub("^label", "", colnames(design))
  fit <- lmFit(expr_mat, design)

  cont <- makeContrasts(
    A_vs_C = A - Control,
    B_vs_C = B - Control,
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, cont))

  ttA <- topTable(fit2, coef = "A_vs_C", number = Inf, sort.by = "none")
  ttB <- topTable(fit2, coef = "B_vs_C", number = Inf, sort.by = "none")

  data.table(
    gene = rownames(ttA),
    delta_1 = ttA$logFC,
    delta_2 = ttB$logFC
  )
}

run_pooled_limma_timecourse <- function(expr_mat, sample_info) {
  design <- model.matrix(~ 0 + group, data = sample_info)
  colnames(design) <- sub("^group", "", colnames(design))
  fit <- lmFit(expr_mat, design)
  cont <- makeContrasts(Treatment - Control, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), pooled_logFC = tt$logFC)
}

run_timeaware_limma_timecourse <- function(expr_mat, sample_info) {
  design <- model.matrix(~ 0 + label, data = sample_info)
  colnames(design) <- sub("^label", "", colnames(design))
  fit <- lmFit(expr_mat, design)

  cont <- makeContrasts(
    t1_vs_C = t1 - Control,
    t2_vs_C = t2 - Control,
    t3_vs_C = t3 - Control,
    t4_vs_C = t4 - Control,
    t5_vs_C = t5 - Control,
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, cont))

  out_list <- vector("list", length(timepoints))
  coef_names <- c("t1_vs_C", "t2_vs_C", "t3_vs_C", "t4_vs_C", "t5_vs_C")
  tp_names   <- c("t1", "t2", "t3", "t4", "t5")

  for (i in seq_along(coef_names)) {
    tt <- topTable(fit2, coef = coef_names[i], number = Inf, sort.by = "none")
    out_list[[i]] <- data.table(
      gene = rownames(tt),
      timepoint = tp_names[i],
      logFC = tt$logFC
    )
  }

  rbindlist(out_list)
}

# =========================================================
# Synthetic generators
# =========================================================
simulate_subgroup_dataset <- function(noise_sd) {
  control_cols <- paste0("C", seq_len(n_rep))
  A_cols <- paste0("A", seq_len(n_rep))
  B_cols <- paste0("B", seq_len(n_rep))
  all_cols <- c(control_cols, A_cols, B_cols)

  expr <- matrix(NA_real_, nrow = n_genes, ncol = length(all_cols),
                 dimnames = list(genes, all_cols))

  expr[, control_cols] <- replicate(n_rep, rnorm(n_genes, 0, noise_sd))

  idx_null     <- truth$scenario == "null"
  idx_coherent <- truth$scenario == "coherent"
  idx_opposing <- truth$scenario == "opposing"
  idx_phase    <- truth$scenario == "phase"

  A_mean <- numeric(n_genes)
  A_mean[idx_null]     <- 0
  A_mean[idx_coherent] <- 1.0
  A_mean[idx_opposing] <- 1.5
  A_mean[idx_phase]    <- 0

  B_mean <- numeric(n_genes)
  B_mean[idx_null]     <- 0
  B_mean[idx_coherent] <- 1.2
  B_mean[idx_opposing] <- -1.5
  B_mean[idx_phase]    <- 0

  expr[, A_cols] <- replicate(n_rep, rnorm(n_genes, A_mean, noise_sd))
  expr[, B_cols] <- replicate(n_rep, rnorm(n_genes, B_mean, noise_sd))

  sample_info <- data.table(
    sample = all_cols,
    group = c(rep("Control", n_rep), rep("Treatment", 2 * n_rep)),
    label = c(rep("Control", n_rep), rep("A", n_rep), rep("B", n_rep))
  )

  list(expr = expr, sample_info = sample_info)
}

simulate_timecourse_dataset <- function(noise_sd) {
  control_cols <- paste0("C", seq_len(n_rep))
  tp_cols <- unlist(lapply(timepoints, function(tp) paste0(tp, "_", seq_len(n_rep))))
  all_cols <- c(control_cols, tp_cols)

  expr <- matrix(NA_real_, nrow = n_genes, ncol = length(all_cols),
                 dimnames = list(genes, all_cols))

  expr[, control_cols] <- replicate(n_rep, rnorm(n_genes, 0, noise_sd))

  idx_null     <- truth$scenario == "null"
  idx_coherent <- truth$scenario == "coherent"
  idx_opposing <- truth$scenario == "opposing"
  idx_phase    <- truth$scenario == "phase"

  phase_means <- c(t1 = 1.2, t2 = 0.8, t3 = 0.0, t4 = -0.8, t5 = -1.1)

  for (tp in timepoints) {
    cols <- paste0(tp, "_", seq_len(n_rep))
    tp_mean <- numeric(n_genes)

    tp_mean[idx_null]     <- 0
    tp_mean[idx_coherent] <- 1.0
    tp_mean[idx_opposing] <- 0
    tp_mean[idx_phase]    <- phase_means[[tp]]

    expr[, cols] <- replicate(n_rep, rnorm(n_genes, tp_mean, noise_sd))
  }

  sample_info <- data.table(
    sample = all_cols,
    group = c(rep("Control", n_rep), rep("Treatment", length(tp_cols))),
    label = c(rep("Control", n_rep), rep(rep(timepoints, each = n_rep), times = 1))
  )

  list(expr = expr, sample_info = sample_info)
}

# =========================================================
# Repeated simulations
# =========================================================
subgroup_runs <- vector("list", n_sim)
timecourse_runs <- vector("list", n_sim)

for (sim in seq_len(n_sim)) {

  # subgroup
  sg <- simulate_subgroup_dataset(noise_sd)
  pooled_sg <- run_pooled_limma_subgroup(sg$expr, sg$sample_info)
  strat_sg  <- run_stratified_limma_subgroup(sg$expr, sg$sample_info)

  sg_sum <- merge(pooled_sg, strat_sg, by = "gene")
  sg_sum <- merge(sg_sum, truth, by = "gene")

  sg_summary <- sg_sum[
    scenario %in% c("null", "coherent", "opposing"),
    .(
      pooled_abs_logFC     = mean(abs(pooled_logFC), na.rm = TRUE),
      structured_abs_logFC = mean((abs(delta_1) + abs(delta_2)) / 2, na.rm = TRUE)
    ),
    by = scenario
  ]

  sg_summary[, `:=`(sim = sim, panel = "Subgroup")]
  subgroup_runs[[sim]] <- sg_summary

  # timecourse
  tc <- simulate_timecourse_dataset(noise_sd)
  pooled_tc <- run_pooled_limma_timecourse(tc$expr, tc$sample_info)
  aware_tc  <- run_timeaware_limma_timecourse(tc$expr, tc$sample_info)

  aware_wide <- dcast(aware_tc, gene ~ timepoint, value.var = "logFC")
  tc_sum <- merge(pooled_tc, aware_wide, by = "gene")
  tc_sum <- merge(tc_sum, truth, by = "gene")

  tc_summary <- tc_sum[
    scenario %in% c("null", "coherent", "phase"),
    .(
      pooled_abs_logFC     = mean(abs(pooled_logFC), na.rm = TRUE),
      structured_abs_logFC = mean(rowMeans(abs(as.matrix(.SD))), na.rm = TRUE)
    ),
    .SDcols = timepoints,
    by = scenario
  ]

  tc_summary[, `:=`(sim = sim, panel = "Timecourse")]
  timecourse_runs[[sim]] <- tc_summary
}

subgroup_dt   <- rbindlist(subgroup_runs, use.names = TRUE, fill = TRUE)
timecourse_dt <- rbindlist(timecourse_runs, use.names = TRUE, fill = TRUE)

# =========================================================
# Prepare plotting table
# =========================================================
sg_plot <- melt(
  subgroup_dt,
  id.vars = c("sim", "panel", "scenario"),
  measure.vars = c("pooled_abs_logFC", "structured_abs_logFC"),
  variable.name = "method",
  value.name = "value"
)

tc_plot <- melt(
  timecourse_dt,
  id.vars = c("sim", "panel", "scenario"),
  measure.vars = c("pooled_abs_logFC", "structured_abs_logFC"),
  variable.name = "method",
  value.name = "value"
)

plot_dt <- rbindlist(list(sg_plot, tc_plot), use.names = TRUE)

plot_dt[, method := fifelse(method == "pooled_abs_logFC", "Pooled", "Structured")]
plot_dt[, panel := factor(panel, levels = c("Subgroup", "Timecourse"))]

# ---- panel-specific plotting variable ----
plot_dt[, scenario_plot := NA_character_]

plot_dt[panel == "Subgroup" & scenario == "null",     scenario_plot := "sub_null"]
plot_dt[panel == "Subgroup" & scenario == "coherent", scenario_plot := "sub_coherent"]
plot_dt[panel == "Subgroup" & scenario == "opposing", scenario_plot := "sub_opposing"]

plot_dt[panel == "Timecourse" & scenario == "null",     scenario_plot := "time_null"]
plot_dt[panel == "Timecourse" & scenario == "coherent", scenario_plot := "time_coherent"]
plot_dt[panel == "Timecourse" & scenario == "phase",    scenario_plot := "time_phase"]

plot_dt[, scenario_plot := factor(
  scenario_plot,
  levels = c("sub_null", "sub_coherent", "sub_opposing",
             "time_null", "time_coherent", "time_phase")
)]

summary_dt <- plot_dt[, .(
  mean_value = mean(value, na.rm = TRUE),
  sd_value   = sd(value, na.rm = TRUE)
), by = .(panel, scenario_plot, method)]

summary_dt[, panel := factor(panel, levels = c("Subgroup", "Timecourse"))]

# =========================================================
# Plot
# =========================================================
p <- ggplot(summary_dt, aes(x = scenario_plot, y = mean_value, fill = method)) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.72,
    color = "black",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
    position = position_dodge(width = 0.78),
    width = 0.16,
    linewidth = 0.55
  ) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    values = c("Pooled" = "#1f77b4", "Structured" = "#ff7f0e")
  ) +
  scale_x_discrete(
  labels = c(
    sub_null = "null",
    sub_coherent = "coherent",
    sub_opposing = "opposing",
    time_null = "null",
    time_coherent = "coherent",
    time_phase = "phase"
  ),
  drop = TRUE
) +
  labs(
    title = "Pooled analysis attenuates signal under structured heterogeneity",
    x = NULL,
    y = "Mean absolute log2 fold change"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    strip.text = element_text(size = 24),
    legend.title = element_blank(),
    legend.text = element_text(size = 22),
    legend.position = "top",
    panel.grid = element_blank(),
    plot.margin = margin(12, 20, 12, 12)
  ) +
  coord_cartesian(
    ylim = c(0, max(summary_dt$mean_value + summary_dt$sd_value, na.rm = TRUE) + 0.15)
  )

try(
  ggsave(
    filename = file.path(out_dir, "figure2B_final_lock.png"),
    plot = p,
    width = 16,
    height = 10,
    dpi = 300,
    device = "png"
  ),
  silent = TRUE
)

ggsave(
  filename = file.path(out_dir, "figure2B_final_lock.pdf"),
  plot = p,
  width = 16,
  height = 10,
  device = "pdf"
)

if (!file.exists(file.path(out_dir, "figure2B_final_lock.pdf"))) {
  stop("PDF figure export failed: ", file.path(out_dir, "figure2B_final_lock.pdf"))
}

# =========================================================
# Save summary table
# =========================================================
fwrite(summary_dt, file.path(out_dir, "figure2B_summary_stats.tsv"), sep = "\t")

cat("Figure 2B locked output written to:\n")
cat(file.path(out_dir, "figure2B_final_lock.png"), "\n")
cat(file.path(out_dir, "figure2B_summary_stats.tsv"), "\n")
