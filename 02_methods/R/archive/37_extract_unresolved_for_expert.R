#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
in_file <- file.path(
  nm2_root,
  "04_results", "figure5e_FINAL", "figure5e_raw_annotation_table.tsv"
)

out_dir <- file.path(
  nm2_root,
  "04_results", "figure5e_FINAL", "expert_annotation"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "figure5e_unresolved_for_expert.tsv")

if (!file.exists(in_file)) stop("Missing file: ", in_file)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
cat("Loading annotation table...\n")
ann <- read_tsv(in_file, show_col_types = FALSE)

required_cols <- c("gene", "protein_id", "raw_annotation", "broad_group_figure")
miss <- setdiff(required_cols, names(ann))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

# ------------------------------------------------------------
# Extract unresolved
# ------------------------------------------------------------
unresolved <- ann %>%
  filter(broad_group_figure == "Complex proteins, unresolved class") %>%
  arrange(desc(instability_score)) %>%
  select(
    gene,
    protein_id,
    delta_pool,
    instability_score,
    match_type,
    ncbi_protein_desc,
    pfam_name,
    pfam_acc,
    pfam_desc,
    pfam_hit_count,
    interpro_analysis,
    interpro_sig_acc,
    interpro_sig_desc,
    interpro_acc,
    interpro_desc,
    interpro_hit_count,
    raw_annotation
  ) %>%
  mutate(
    expert_annotation = "",     # for expert to fill
    expert_category = ""        # optional structured label
  )

write_tsv(unresolved, out_file)

cat("\nExtracted unresolved proteins:", nrow(unresolved), "\n")
cat("Saved to:\n", out_file, "\n")
