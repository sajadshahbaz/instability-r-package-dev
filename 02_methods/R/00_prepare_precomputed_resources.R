#!/usr/bin/env Rscript

root <- Sys.getenv("NM2_ROOT")
if (root == "") root <- normalizePath(getwd(), mustWork = TRUE)

copy_required <- function(from, to) {
  if (!file.exists(from)) stop("Missing bundled resource: ", from)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  file.copy(from, to, overwrite = TRUE)
  if (!file.exists(to)) stop("Failed to copy resource to: ", to)
}

copy_required(
  file.path(root, "01_data/precomputed/figure5e_annotation/figure5e_pfam.domtblout"),
  file.path(root, "04_results/figure5e_global_functional_structure/annotation_from_scratch/figure5e_pfam.domtblout")
)

copy_required(
  file.path(root, "01_data/precomputed/figure5e_annotation/figure5e_interproscan.tsv"),
  file.path(root, "04_results/figure5e_global_functional_structure/annotation_from_scratch/figure5e_interproscan.tsv")
)

copy_required(
  file.path(root, "01_data/precomputed/figure6/Figure6_variance_disentanglement_FINAL.pdf"),
  file.path(root, "04_results/figure6_variance_disentanglement_FINAL/Figure6_variance_disentanglement_FINAL.pdf")
)

copy_required(
  file.path(root, "01_data/precomputed/figure6/figure6_source_table.tsv"),
  file.path(root, "04_results/figure6_variance_disentanglement_FINAL/figure6_source_table.tsv")
)

cat("Precomputed manuscript resources prepared.\n")
