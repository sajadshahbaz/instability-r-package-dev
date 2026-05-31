#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT is not set")
base_dir <- file.path(nm2_root, "04_results", "figure5c_ncbi_clean")

gene_map <- file.path(base_dir, "figure5c_genes_matched_to_ncbi.csv")
fasta_file <- file.path(base_dir, "figure5c_matched_ncbi_proteins.faa")
merged_file <- file.path(base_dir, "figure5c_final/figure5c_fresh_annotation_merged.csv")

cat("\n=== STEP 1: Load gene mapping ===\n")

gm <- read_csv(gene_map, show_col_types = FALSE)

cat("Total rows:", nrow(gm), "\n")
cat("Unique genes:", n_distinct(gm$gene), "\n")
cat("Unique protein_ids:", n_distinct(gm$protein_id), "\n")

# Check duplicates
dup_genes <- gm %>% count(gene) %>% filter(n > 1)
dup_prot  <- gm %>% count(protein_id) %>% filter(n > 1)

cat("Duplicated genes:", nrow(dup_genes), "\n")
cat("Duplicated protein_ids:", nrow(dup_prot), "\n")

if (nrow(dup_genes) > 0) {
  cat("Example duplicated genes:\n")
  print(head(dup_genes))
}

if (nrow(dup_prot) > 0) {
  cat("Example duplicated proteins:\n")
  print(head(dup_prot))
}

# Missing values
cat("Missing gene:", sum(is.na(gm$gene)), "\n")
cat("Missing protein_id:", sum(is.na(gm$protein_id)), "\n")

# ---------------------------------------------------
cat("\n=== STEP 2: Check FASTA consistency ===\n")

fasta_lines <- readLines(fasta_file)
headers <- fasta_lines[grepl("^>", fasta_lines)]

extract_id <- function(x) {
  x <- sub("^>", "", x)
  strsplit(x, " ")[[1]][1]
}

fasta_ids <- sapply(headers, extract_id)

cat("FASTA entries:", length(fasta_ids), "\n")
cat("Unique FASTA IDs:", length(unique(fasta_ids)), "\n")

# Missing from FASTA
missing_in_fasta <- setdiff(gm$protein_id, fasta_ids)
extra_in_fasta   <- setdiff(fasta_ids, gm$protein_id)

cat("Missing protein_ids in FASTA:", length(missing_in_fasta), "\n")
cat("Extra FASTA entries not in mapping:", length(extra_in_fasta), "\n")

if (length(missing_in_fasta) > 0) {
  cat("Example missing IDs:\n")
  print(head(missing_in_fasta))
}

# ---------------------------------------------------
cat("\n=== STEP 3: Check merged annotation consistency ===\n")

df <- read_csv(merged_file, show_col_types = FALSE)

cat("Merged rows:", nrow(df), "\n")
cat("Unique protein_ids:", n_distinct(df$protein_id), "\n")

# Check join consistency
missing_in_annotation <- setdiff(gm$protein_id, df$protein_id)

cat("Protein IDs missing in annotation:", length(missing_in_annotation), "\n")

if (length(missing_in_annotation) > 0) {
  print(head(missing_in_annotation))
}

# ---------------------------------------------------
cat("\n=== STEP 4: ID formatting issues ===\n")

cat("Example protein IDs:\n")
print(head(gm$protein_id, 10))

# detect prefixes
has_lcl <- any(grepl("^lcl\\|", gm$protein_id))
cat("Contains lcl| prefix:", has_lcl, "\n")

has_version <- any(grepl("\\.[0-9]+$", gm$protein_id))
cat("Contains version suffix (.1, .2):", has_version, "\n")

cat("\n=== VALIDATION COMPLETE ===\n")
