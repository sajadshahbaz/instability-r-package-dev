#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Figure 5C — Match NM2 genes directly to NCBI reference
# ============================================================

# -----------------------------
# Paths
# -----------------------------
nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT is not set")
base_topn <- file.path(nm2_root, "04_results", "figure5b_topN_overlap")

inst_file <- file.path(base_topn, "instability_only_genes.tsv")
pool_file <- file.path(base_topn, "pooled_only_genes.tsv")
over_file <- file.path(base_topn, "overlap_genes.tsv")

work_dir <- file.path(nm2_root, "04_results", "figure5c_ncbi_clean")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

ref_file <- file.path(work_dir, "figure5c_ncbi_reference_map.csv")

out_master <- file.path(work_dir, "figure5c_gene_sets_master.csv")
out_match <- file.path(work_dir, "figure5c_genes_matched_to_ncbi.csv")
out_audit <- file.path(work_dir, "figure5c_match_audit.csv")
out_unmatched <- file.path(work_dir, "figure5c_unmatched_genes.csv")

# -----------------------------
# Helpers
# -----------------------------
normalize_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_.-]", "", x)
  x <- gsub("(?i)^rvy_", "RvY_", x, perl = TRUE)
  x
}

base_gene <- function(x) {
  x <- normalize_gene(x)
  gsub("-\\d+$", "", x)
}

read_gene_file <- function(fp, label) {
  if (!file.exists(fp)) stop("Missing gene file: ", fp)
  df <- read_tsv(fp, show_col_types = FALSE)
  if (!"gene" %in% colnames(df)) stop("Expected 'gene' column in: ", fp)

  df %>%
    transmute(
      gene = normalize_gene(gene),
      category = label,
      gene_base = base_gene(gene)
    ) %>%
    filter(gene != "") %>%
    distinct()
}

# -----------------------------
# Load NM2 genes
# -----------------------------
inst <- read_gene_file(inst_file, "Instability-only")
pool <- read_gene_file(pool_file, "Pooled-only")
over <- read_gene_file(over_file, "Overlap")

master <- bind_rows(inst, pool, over) %>%
  distinct()

write_csv(master, out_master)

# -----------------------------
# Load reference map
# -----------------------------
if (!file.exists(ref_file)) stop("Missing reference map: ", ref_file)
ref <- read_csv(ref_file, show_col_types = FALSE)

# -----------------------------
# Matching hierarchy
# -----------------------------
m_exact_gene <- master %>%
  left_join(ref, by = c("gene" = "gene")) %>%
  filter(!is.na(protein_id)) %>%
  mutate(match_type = "exact_gene")

matched_genes_1 <- unique(m_exact_gene$gene)

m_exact_locus <- master %>%
  filter(!(gene %in% matched_genes_1)) %>%
  left_join(ref, by = c("gene" = "locus_tag")) %>%
  filter(!is.na(protein_id)) %>%
  mutate(match_type = "exact_locus")

matched_genes_2 <- unique(c(matched_genes_1, m_exact_locus$gene))

m_base_gene <- master %>%
  filter(!(gene %in% matched_genes_2)) %>%
  left_join(ref, by = c("gene_base" = "gene_base")) %>%
  filter(!is.na(protein_id)) %>%
  mutate(match_type = "base_gene")

matched_genes_3 <- unique(c(matched_genes_2, m_base_gene$gene))

m_base_locus <- master %>%
  filter(!(gene %in% matched_genes_3)) %>%
  left_join(ref, by = c("gene_base" = "locus_base")) %>%
  filter(!is.na(protein_id)) %>%
  mutate(match_type = "base_locus")

matched <- bind_rows(
  m_exact_gene,
  m_exact_locus,
  m_base_gene,
  m_base_locus
) %>%
  distinct(gene, category, .keep_all = TRUE)

unmatched <- master %>%
  filter(!(gene %in% matched$gene))

write_csv(matched, out_match)
write_csv(unmatched, out_unmatched)

audit <- matched %>%
  group_by(category) %>%
  summarise(
    matched = n(),
    .groups = "drop"
  ) %>%
  right_join(
    master %>% count(category, name = "total"),
    by = "category"
  ) %>%
  mutate(
    matched = coalesce(matched, 0L),
    unmatched = total - matched,
    matched_pct = round(100 * matched / total, 1)
  )

write_csv(audit, out_audit)

cat("✅ Total input rows:", nrow(master), "\n")
cat("✅ Matched rows:", nrow(matched), "\n")
cat("✅ Unmatched rows:", nrow(unmatched), "\n")
cat("✅ Match type counts:\n")
print(table(matched$match_type, useNA = "ifany"))
print(audit)
