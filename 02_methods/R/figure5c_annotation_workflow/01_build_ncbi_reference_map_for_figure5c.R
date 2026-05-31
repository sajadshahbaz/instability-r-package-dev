#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Figure 5C — Build direct NCBI reference map
# Base-R FASTA parsing version
# ============================================================

# -----------------------------
# Paths
# -----------------------------
nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT is not set")
out_dir <- file.path(nm2_root, "04_results", "figure5c_ncbi_clean")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ncbi_root <- file.path(
  nm2_root,
  "01_data", "raw", "ncbi_dataset", "ncbi_dataset", "data", "GCA_001949185.1"
)

cds_fa  <- file.path(ncbi_root, "cds_from_genomic.fna")
prot_fa <- file.path(ncbi_root, "protein.faa")

meta_file <- file.path(out_dir, "REFERENCE_METADATA.txt")

# -----------------------------
# Helpers
# -----------------------------
extract_field <- function(x, key) {
  m <- str_match(x, paste0("\\[", key, "=([^\\]]+)\\]"))[,2]
  m
}

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

read_fasta_headers <- function(fp) {
  if (!file.exists(fp)) stop("Missing FASTA: ", fp)
  x <- readLines(fp, warn = FALSE)
  x[startsWith(x, ">")]
}

# -----------------------------
# Input checks
# -----------------------------
if (!file.exists(cds_fa))  stop("Missing CDS FASTA: ", cds_fa)
if (!file.exists(prot_fa)) stop("Missing protein FASTA: ", prot_fa)

# -----------------------------
# Read CDS headers
# -----------------------------
cds_hdr_raw <- read_fasta_headers(cds_fa)
cds_hdr <- sub("^>", "", cds_hdr_raw)

cds_map <- tibble(
  cds_header = cds_hdr,
  cds_seqname = sub("\\s.*$", "", cds_hdr),
  gene = normalize_gene(extract_field(cds_hdr, "gene")),
  locus_tag = normalize_gene(extract_field(cds_hdr, "locus_tag")),
  protein_desc_from_cds = extract_field(cds_hdr, "protein"),
  protein_id = extract_field(cds_hdr, "protein_id"),
  gene_base = base_gene(extract_field(cds_hdr, "gene")),
  locus_base = base_gene(extract_field(cds_hdr, "locus_tag"))
)

# -----------------------------
# Read protein headers
# -----------------------------
prot_hdr_raw <- read_fasta_headers(prot_fa)
prot_hdr <- sub("^>", "", prot_hdr_raw)

prot_desc <- sub("^[^ ]+\\s*", "", prot_hdr)
prot_desc <- str_replace(prot_desc, "\\s*\\[Ramazzottius varieornatus\\]$", "")

prot_map <- tibble(
  protein_header = prot_hdr,
  protein_id = sub("\\s.*$", "", prot_hdr),
  protein_desc_from_protein = prot_desc
)

# -----------------------------
# Join
# -----------------------------
ref_map <- cds_map %>%
  left_join(prot_map, by = "protein_id") %>%
  mutate(
    protein_present_in_faa = !is.na(protein_header)
  ) %>%
  distinct()

write_csv(ref_map, file.path(out_dir, "figure5c_ncbi_reference_map.csv"))

# -----------------------------
# Metadata
# -----------------------------
meta_lines <- c(
  "Figure 5C NCBI reference metadata",
  "Organism: Ramazzottius varieornatus",
  "Source: NCBI Datasets taxonomy page https://www.ncbi.nlm.nih.gov/datasets/taxonomy/947166/",
  paste0("NCBI root: ", ncbi_root),
  paste0("CDS FASTA: ", cds_fa),
  paste0("Protein FASTA: ", prot_fa),
  paste0("Build date: ", Sys.time())
)
writeLines(meta_lines, meta_file)

cat("✅ CDS entries:", nrow(cds_map), "\n")
cat("✅ Protein entries:", nrow(prot_map), "\n")
cat("✅ Reference map rows:", nrow(ref_map), "\n")
cat("✅ Rows with matched protein accession in protein.faa:", sum(ref_map$protein_present_in_faa), "\n")
