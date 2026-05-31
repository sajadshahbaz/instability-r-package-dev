# NM2 instability-signal master runner
# Reproduces the main manuscript workflow from organized repository scripts.

message("Starting NM2 instability-signal pipeline...")

root <- normalizePath(getwd(), mustWork = TRUE)
Sys.setenv(NM2_ROOT = root)
message("Repository root: ", root)
message("NM2_ROOT set to: ", Sys.getenv("NM2_ROOT"))
run_script <- function(path) {
  full_path <- file.path(root, path)

  if (!file.exists(full_path)) {
    stop("Missing script: ", full_path)
  }

  message("\n--- Running: ", path, " ---")

  res <- system2(
    "Rscript",
    c("--vanilla", full_path),
    stdout = "",
    stderr = ""
  )

  if (!identical(res, 0L)) {
    stop("Script failed: ", path)
  }
}

# Core pipeline
core_scripts <- c(
  "02_methods/R/00_prepare_precomputed_resources.R",
  "02_methods/R/core_pipeline/00_io.R",
  "02_methods/R/core_pipeline/01_make_homogeneous_null_controls.R",
  "02_methods/R/core_pipeline/01_pool_vs_strat_rebuild_uv_des.R",
  "02_methods/R/core_pipeline/02_plots_pool_vs_strat.R",
  "02_methods/R/core_pipeline/03_tail_divergence.R",
  "02_methods/R/core_pipeline/04_make_core_table.R",
  "02_methods/R/core_pipeline/05_gamma_instability_from_timepoints.R",
  "02_methods/R/core_pipeline/07_gate3_instability_vs_variance.R",
  "02_methods/R/core_pipeline/08_gate4_regime_conflict_enrichment.R"
  )

# Figure generation
figure_scripts <- c(
  "02_methods/R/figure_generation/figure2_pooled_failure.R",
  "02_methods/R/figure_generation/figure2_instability_enrichment.R",
  "02_methods/R/figure_generation/figure3_uv_quadrant.R",
  "02_methods/R/figure_generation/figure3_instability_vs_variance.R",
  "02_methods/R/figure_generation/figure3_instability_vs_disagreement.R",
  "02_methods/R/figure_generation/figure3_gene_examples.R",
  "02_methods/R/figure_generation/figure4_des_support.R",
  "02_methods/R/figure_generation/figure4_gam_temporal_extension.R",
  "02_methods/R/figure_generation/figure4_cross_condition_panel.R",
  "02_methods/R/figure_generation/figure5_rank_displacement.R",
  "02_methods/R/figure_generation/figure5_topN_overlap.R",
  "02_methods/R/figure_generation/figure5c_selection_stability.R",
  "02_methods/R/figure5c_annotation_workflow/01_build_ncbi_reference_map_for_figure5c.R",
  "02_methods/R/figure5c_annotation_workflow/02_match_figure5c_genes_to_ncbi_reference.R",
  "02_methods/R/figure5c_annotation_workflow/03_extract_matched_ncbi_proteins_for_figure5c.R",
  "02_methods/R/figure5c_annotation_workflow/04_merge_annotation_and_make_figure5c.R",
  "02_methods/R/figure5c_annotation_workflow/00_validate_mapping_integrity.R",
  "02_methods/R/figure_generation/figure5d_curated_examples_table.R",
  "02_methods/R/figure_generation/figure5e_functional_composition.R"
)

all_scripts <- c(core_scripts, figure_scripts)

for (script in all_scripts) {
  run_script(script)
}

message("\nNM2 pipeline completed successfully.")
