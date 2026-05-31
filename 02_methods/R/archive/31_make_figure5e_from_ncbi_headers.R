#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
uv_file <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")

ncbi_dir <- "/media/saji/5E06441D0643F5152/nature_method2/NM2_instability_signal/01_data/raw/ncbi_dataset/ncbi_dataset/data/GCA_001949185.1"
cds_fna  <- file.path(ncbi_dir, "cds_from_genomic.fna")
prot_faa <- file.path(ncbi_dir, "protein.faa")

out_dir <- file.path(nm2_root, "04_results", "figure5e_ncbi_header_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

masked_out   <- file.path(out_dir, "figure5e_masked_genes.tsv")
map_out      <- file.path(out_dir, "figure5e_gene_to_protein_from_ncbi.tsv")
counts_out   <- file.path(out_dir, "figure5e_category_counts.tsv")
annot_out    <- file.path(out_dir, "figure5e_ncbi_annotation_table.tsv")
coverage_out <- file.path(out_dir, "figure5e_annotation_coverage.tsv")
plot_png     <- file.path(out_dir, "figure5e_ncbi_header_annotation.png")
plot_pdf     <- file.path(out_dir, "figure5e_ncbi_header_annotation.pdf")

for (fp in c(uv_file, cds_fna, prot_faa)) {
  if (!file.exists(fp)) stop("Missing file: ", fp)
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
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

clean_txt <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

extract_tag <- function(header, tag) {
  m <- str_match(header, paste0("\\[", tag, "=([^\\]]+)\\]"))
  m[, 2]
}

assign_category <- function(desc) {
  desc2 <- tolower(clean_txt(desc))
  desc2 <- gsub("\\s+", " ", desc2)

  if (
    desc2 == "" ||
    grepl("^hypothetical protein( .*)?$", desc2) ||
    grepl("^uncharacterized protein( .*)?$", desc2) ||
    grepl("^predicted protein( .*)?$", desc2)
  ) {
    return("Uncharacterized/Hypothetical")
  }

  if (grepl("gpcr|receptor|kinase|phosphatase|g-protein|signal", desc2)) {
    return("Signaling/Receptor")
  }
  if (grepl("protease|peptidase|trypsin|serine protease|cathepsin|ubiquitin", desc2)) {
    return("Proteolysis")
  }
  if (grepl("transport|channel|carrier|transporter|membrane", desc2)) {
    return("Transport/Membrane")
  }
  if (grepl("transcription|dna-binding|rna-binding|regulator|chromatin", desc2)) {
    return("Regulation/Nucleic-acid associated")
  }
  if (grepl("dehydrogenase|synthetase|hydrolase|transferase|oxidoreductase|enzyme|metabol", desc2)) {
    return("Metabolism/Enzyme")
  }
  if (grepl("actin|tubulin|cytoskeleton|structural|motor protein|calycin", desc2)) {
    return("Structural/Other")
  }

  return("Other annotated")
}
# ------------------------------------------------------------
# 1) Formal masked genes from UV full table
# ------------------------------------------------------------
cat("Loading UV full table...\n")
uv <- read_tsv(uv_file, show_col_types = FALSE)

if ("feature" %in% names(uv) && !("gene" %in% names(uv))) {
  uv <- uv %>% rename(gene = feature)
}

req_cols <- c("gene", "delta_pool", "instability_score")
if (!all(req_cols %in% names(uv))) {
  stop("UV table missing required columns: ", paste(setdiff(req_cols, names(uv)), collapse = ", "))
}

uv <- uv %>%
  mutate(
    gene = normalize_gene(gene),
    gene_base = base_gene(gene)
  ) %>%
  filter(gene != "")

q25_pool <- quantile(abs(uv$delta_pool), 0.25, na.rm = TRUE)
q90_inst <- quantile(uv$instability_score, 0.90, na.rm = TRUE)

masked <- uv %>%
  filter(abs(delta_pool) <= q25_pool,
         instability_score >= q90_inst) %>%
  distinct(gene, .keep_all = TRUE)

write_tsv(masked, masked_out)
cat("Masked genes:", nrow(masked), "\n")

# ------------------------------------------------------------
# 2) Parse CDS headers for gene -> protein mapping
# ------------------------------------------------------------
cat("Parsing cds_from_genomic.fna headers...\n")
cds_lines <- readLines(cds_fna, warn = FALSE)
cds_headers <- cds_lines[startsWith(cds_lines, ">")]

cds_map <- tibble(
  raw_header = cds_headers,
  gene = extract_tag(cds_headers, "gene"),
  locus_tag = extract_tag(cds_headers, "locus_tag"),
  protein_desc = extract_tag(cds_headers, "protein"),
  protein_id = extract_tag(cds_headers, "protein_id")
) %>%
  mutate(
    gene = normalize_gene(gene),
    locus_tag = normalize_gene(locus_tag),
    gene_base = base_gene(gene),
    locus_base = base_gene(locus_tag)
  ) %>%
  filter(!is.na(protein_id), protein_id != "") %>%
  distinct(gene, .keep_all = TRUE)

# ------------------------------------------------------------
# 3) Parse protein.faa headers for protein description
# ------------------------------------------------------------
cat("Parsing protein.faa headers...\n")
prot_lines <- readLines(prot_faa, warn = FALSE)
prot_headers <- prot_lines[startsWith(prot_lines, ">")]

parse_protein_header <- function(h) {
  hh <- sub("^>", "", h)
  protein_id <- strsplit(hh, " ")[[1]][1]
  rest <- sub(paste0("^", protein_id, " "), "", hh)

  # remove species tail
  rest2 <- sub("\\s*\\[[^\\]]+\\]\\s*$", "", rest)

  tibble(
    protein_id = protein_id,
    protein_header_desc = clean_txt(rest2)
  )
}

prot_map <- bind_rows(lapply(prot_headers, parse_protein_header)) %>%
  distinct(protein_id, .keep_all = TRUE)

# ------------------------------------------------------------
# 4) Match masked genes to NCBI mapping
# ------------------------------------------------------------
cat("Matching masked genes to NCBI gene/protein annotations...\n")

m1 <- masked %>%
  left_join(cds_map, by = c("gene" = "gene")) %>%
  mutate(match_type = case_when(
    !is.na(protein_id) ~ "exact_gene",
    TRUE ~ NA_character_
  ))

matched1 <- m1 %>% filter(!is.na(protein_id))
matched_genes1 <- unique(matched1$gene)

m2 <- masked %>%
  filter(!(gene %in% matched_genes1)) %>%
  left_join(cds_map, by = c("gene_base" = "gene_base")) %>%
  mutate(match_type = case_when(
    !is.na(protein_id) ~ "base_gene",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(protein_id))

matched <- bind_rows(matched1, m2) %>%
  distinct(gene, .keep_all = TRUE) %>%
  left_join(prot_map, by = "protein_id") %>%
  mutate(
    final_desc = dplyr::case_when(
      !is.na(protein_header_desc) & protein_header_desc != "" ~ protein_header_desc,
      !is.na(protein_desc) & protein_desc != "" ~ protein_desc,
      TRUE ~ ""
    ),
    category = vapply(final_desc, assign_category, character(1))
  )

unmatched <- masked %>%
  filter(!(gene %in% matched$gene)) %>%
  mutate(
    protein_id = NA_character_,
    protein_desc = "",
    protein_header_desc = "",
    final_desc = "",
    category = "Uncharacterized/Hypothetical",
    match_type = "unmatched"
  )

annot_tbl <- bind_rows(
  matched %>%
    select(gene, gene_base, protein_id, final_desc, category, match_type),
  unmatched %>%
    select(gene, gene_base, protein_id, final_desc, category, match_type)
)

write_tsv(annot_tbl, annot_out)
write_tsv(annot_tbl, map_out)

# ------------------------------------------------------------
# 5) Coverage
# ------------------------------------------------------------
coverage <- annot_tbl %>%
  mutate(
    desc2 = tolower(clean_txt(final_desc)),
    desc2 = gsub("\\s+", " ", desc2),
    has_named_annotation =
      desc2 != "" &
      !grepl("^hypothetical protein( .*)?$", desc2) &
      !grepl("^uncharacterized protein( .*)?$", desc2) &
      !grepl("^predicted protein( .*)?$", desc2)
  ) %>%
  summarise(
    total = n(),
    matched_to_protein = sum(!is.na(protein_id)),
    unmatched = sum(is.na(protein_id)),
    named_annotation = sum(has_named_annotation),
    hypothetical_or_blank = total - named_annotation,
    named_annotation_pct = round(100 * named_annotation / total, 1)
  )

write_tsv(coverage, coverage_out)

cat("\nAnnotation coverage summary:\n")
print(coverage)

# ------------------------------------------------------------
# 6) Category counts
# ------------------------------------------------------------
counts <- annot_tbl %>%
  count(category, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

write_tsv(counts, counts_out)

cat("\nCategory counts:\n")
print(counts)

# ------------------------------------------------------------
# 7) Plot
# ------------------------------------------------------------
p <- ggplot(counts, aes(x = reorder(category, n), y = n)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Functional structure of instability-recovered masked genes",
    x = "Broad functional class",
    y = "Number of genes"
  )

ggsave(plot_png, p, width = 7, height = 5)
ggsave(plot_pdf, p, width = 7, height = 5)

cat("\n✅ Figure 5E outputs written to:\n", out_dir, "\n")
