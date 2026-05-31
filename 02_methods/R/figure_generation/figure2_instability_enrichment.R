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

out_dir <- file.path(nm2_root, "04_results", "figure2c")
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
# Helper functions
# =========================================================
panel_flag_fraction <- function(dt) {
  # panel-specific thresholds within one simulation and one structure
  x_cut <- quantile(abs(dt$delta_pool), 0.50, na.rm = TRUE)
  y_cut <- quantile(dt$instability_score, 0.90, na.rm = TRUE)

  dt[, masked := abs(delta_pool) <= x_cut & instability_score >= y_cut]

  dt[, .(
    fraction_masked = mean(masked, na.rm = TRUE)
  ), by = scenario]
}

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
sim_results <- vector("list", n_sim)

for (sim in seq_len(n_sim)) {

  # ---------------- subgroup ----------------
  sg <- simulate_subgroup_dataset(noise_sd)
  pooled_sg <- run_pooled_limma_subgroup(sg$expr, sg$sample_info)
  strat_sg  <- run_stratified_limma_subgroup(sg$expr, sg$sample_info)

  sg_sum <- merge(pooled_sg, strat_sg, by = "gene")
  sg_sum <- merge(sg_sum, truth, by = "gene")

  sg_sum[, `:=`(
    delta_pool = pooled_logFC,
    instability_score = abs(delta_1 - delta_2) *
      (1 + fifelse(sign(delta_1) != sign(delta_2) & delta_1 != 0 & delta_2 != 0, 1, 0))
  )]

  sg_frac <- panel_flag_fraction(copy(sg_sum))
  sg_frac[, source := "subgroup"]

  # ---------------- timecourse ----------------
  tc <- simulate_timecourse_dataset(noise_sd)
  pooled_tc <- run_pooled_limma_timecourse(tc$expr, tc$sample_info)
  aware_tc  <- run_timeaware_limma_timecourse(tc$expr, tc$sample_info)

  aware_wide <- dcast(aware_tc, gene ~ timepoint, value.var = "logFC")
  tc_sum <- merge(pooled_tc, aware_wide, by = "gene")
  tc_sum <- merge(tc_sum, truth, by = "gene")

  tc_sum[, delta_pool := pooled_logFC]
  tc_sum[, instability_score := apply(.SD, 1, function(x) {
    gap <- max(x) - min(x)
    disc <- mean(sign(x) != sign(mean(x)))
    abs(gap) * (1 + disc)
  }), .SDcols = timepoints]

  tc_frac <- panel_flag_fraction(copy(tc_sum))
  tc_frac[, source := "timecourse"]

  # ---------------- combine within simulation ----------------
  # coherent and null exist in both structures; average them within simulation
  # opposing exists only in subgroup; phase only in timecourse
  sim_dt <- rbindlist(list(sg_frac, tc_frac), use.names = TRUE, fill = TRUE)
  sim_dt <- sim_dt[, .(
    fraction_masked = mean(fraction_masked, na.rm = TRUE)
  ), by = scenario]
  sim_dt[, sim := sim]

  sim_results[[sim]] <- sim_dt
}

all_frac <- rbindlist(sim_results, use.names = TRUE, fill = TRUE)

all_frac[, scenario := factor(
  scenario,
  levels = c("null", "coherent", "opposing", "phase")
)]

summary_dt <- all_frac[, .(
  mean_fraction = mean(fraction_masked, na.rm = TRUE),
  sd_fraction   = sd(fraction_masked, na.rm = TRUE)
), by = scenario]

summary_dt[, scenario := factor(
  scenario,
  levels = c("null", "coherent", "opposing", "phase")
)]

# =========================================================
# Plot
# =========================================================
p <- ggplot(summary_dt, aes(x = scenario, y = mean_fraction)) +
  geom_col(
    width = 0.78,
    fill = "grey35",
    color = "black",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(ymin = mean_fraction - sd_fraction, ymax = mean_fraction + sd_fraction),
    width = 0.14,
    linewidth = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf("%.3f", mean_fraction),
      y = mean_fraction + sd_fraction + 0.02
    ),
    size = 6,
    fontface = "bold"
  ) +
  labs(
    title = "Instability highlights regimes where pooled analysis fails",
    x = NULL,
    y = "Fraction of genes with low-to-moderate pooled effect\nand high instability"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", size = 28, hjust = 0.5),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid = element_blank(),
    plot.margin = margin(12, 20, 12, 12)
  ) +
  coord_cartesian(
    ylim = c(
      0,
      max(summary_dt$mean_fraction + summary_dt$sd_fraction, na.rm = TRUE) + 0.15
    )
  )
try(
try(
    ggsave(
      filename = file.path(out_dir, "figure2C_final_lock.png"),
      plot = p,
      width = 12,
      height = 9,
      dpi = 300,
      bg = "white"
    )
,
  silent = TRUE
)
,
  silent = TRUE
)

ggsave(
  filename = file.path(out_dir, "figure2C_final_lock.pdf"),
  plot = p,
  width = 12,
  height = 9,
  bg = "white"
)

# =========================================================
# Save summary table
# =========================================================
fwrite(summary_dt, file.path(out_dir, "figure2C_summary_stats.tsv"), sep = "\t")

cat("Figure 2C locked output written to:\n")
cat(file.path(out_dir, "figure2C_final_lock.png"), "\n")
cat(file.path(out_dir, "figure2C_summary_stats.tsv"), "\n")
