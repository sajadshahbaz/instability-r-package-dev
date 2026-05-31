#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Figure 5C — Extract matched NCBI proteins
# Base-R FASTA parsing version
# ============================================================

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT is not set")
work_dir <- file.path(nm2_root, "04_results", "figure5c_ncbi_clean")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

ncbi_root <- file.path(
  nm2_root,
  "01_data", "raw", "ncbi_dataset", "ncbi_dataset", "data", "GCA_001949185.1"
)
prot_fa <- file.path(ncbi_root, "protein.faa")

match_file <- file.path(work_dir, "figure5c_genes_matched_to_ncbi.csv")
out_fa <- file.path(work_dir, "figure5c_matched_ncbi_proteins.faa")
out_audit <- file.path(work_dir, "figure5c_protein_extraction_audit.csv")

if (!file.exists(match_file)) stop("Missing matched gene file: ", match_file)
if (!file.exists(prot_fa)) stop("Missing protein FASTA: ", prot_fa)

matched <- read_csv(match_file, show_col_types = FALSE) %>%
  filter(!is.na(protein_id)) %>%
  distinct(gene, category, protein_id)

target_ids <- unique(matched$protein_id)

# -----------------------------
# FASTA parser
# -----------------------------
x <- readLines(prot_fa, warn = FALSE)
hdr_idx <- which(startsWith(x, ">"))

if (length(hdr_idx) == 0) stop("No FASTA headers found in: ", prot_fa)

end_idx <- c(hdr_idx[-1] - 1, length(x))

headers <- sub("^>", "", x[hdr_idx])
prot_ids <- sub("\\s.*$", "", headers)

keep <- prot_ids %in% target_ids

# Ensure output target is writable and not an existing directory
if (dir.exists(out_fa)) {
  stop("Output FASTA path is a directory, not a file: ", out_fa)
}
dir.create(dirname(out_fa), recursive = TRUE, showWarnings = FALSE)

if (file.exists(out_fa)) {
  file.remove(out_fa)
}

written_ids <- character(0)
kept_idx <- which(keep)

if (length(kept_idx) == 0) {
  stop("No target proteins found in protein FASTA. Check protein_id matching.")
}

con <- file(out_fa, open = "wt")
if (!isOpen(con, "write")) {
  stop("Could not open output FASTA for writing: ", out_fa)
}
on.exit({
  try(close(con), silent = TRUE)
}, add = TRUE)

for (i in kept_idx) {
  block <- x[hdr_idx[i]:end_idx[i]]
  writeLines(block, con = con, sep = "\n")
  written_ids <- c(written_ids, prot_ids[i])
}

flush(con)
close(con)

audit <- matched %>%
  mutate(protein_extracted = protein_id %in% written_ids)

write_csv(audit, out_audit)

cat("✅ Target protein IDs:", length(target_ids), "\n")
cat("✅ Extracted protein FASTA entries:", length(unique(written_ids)), "\n")
cat("✅ Audit file written:", out_audit, "\n")
cat("✅ FASTA written:", out_fa, "\n")
