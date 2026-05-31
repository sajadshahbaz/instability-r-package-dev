#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(data.table)
})

# ============================================================
# Figure 5E — Global functional structure of masked genes
# Uses:
#   - UV full gene-level table
#   - existing Figure 5C old annotation resources
# Outputs:
#   - masked gene table
#   - matched annotation table
#   - category counts
#   - barplot
# ============================================================

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT is not set")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
uv_file <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")

work_dir <- file.path(nm2_root, "04_results", "figure5e_global_functional_structure")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# Reuse old Figure 5C resources
fig5c_ref_dir <- file.path(nm2_root, "04_results", "old", "figure5c_ncbi_clean")
ref_file      <- file.path(fig5c_ref_dir, "figure5c_ncbi_reference_map.csv")
annot_file    <- file.path(fig5c_ref_dir, "figure5c_final", "figure5c_fresh_annotation_merged.csv")

masked_out       <- file.path(work_dir, "figure5e_masked_genes.tsv")
matched_out      <- file.path(work_dir, "figure5e_masked_genes_matched.csv")
annot_merged_out <- file.path(work_dir, "figure5e_annotation_merged.csv")
counts_out       <- file.path(work_dir, "figure5e_category_counts.tsv")
coverage_out     <- file.path(work_dir, "figure5e_annotation_coverage.tsv")
plot_png         <- file.path(work_dir, "figure5e_global_functional_structure.png")
plot_pdf         <- file.path(work_dir, "figure5e_global_functional_structure.pdf")
unmatched_out    <- file.path(work_dir, "figure5e_unmatched_genes.tsv")

for (fp in c(uv_file, ref_file, annot_file)) {
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
  x[x %in% c("-", "/", "NA")] <- ""
  trimws(x)
}

collapse_unique <- function(x) {
  x <- clean_txt(x)
  x <- x[x != ""]
  if (!length(x)) return("")
  paste(unique(x), collapse = "; ")
}

# Broad category assignment
assign_category <- function(txt) {
  txt <- clean_txt(txt)
  txt <- tolower(txt)

  if (txt == "") return("Uncharacterized")

  if (grepl("gpcr|receptor|signal transduction|kinase|phosphatase|g-protein|small gtpase|arf|sar", txt)) {
    return("Signaling/Receptor")
  }
  if (grepl("transcription|transcription factor|dna-binding|rna-binding|regulator|chromatin", txt)) {
    return("Regulation")
  }
  if (grepl("transport|channel|carrier|transporter|abc transporter|membrane", txt)) {
    return("Transport/Membrane")
  }
  if (grepl("protease|peptidase|ubiquitin|proteolysis|cathepsin|trypsin|serine protease", txt)) {
    return("Proteolysis")
  }
  if (grepl("heat shock|oxidative|stress|desiccation|response to stress", txt)) {
    return("Stress response")
  }
  if (grepl("dna repair|repair|recombination|replication", txt)) {
    return("DNA repair / genome maintenance")
  }
  if (grepl("cytoskeleton|actin|tubulin|structural|motor protein|calycin", txt)) {
    return("Structural / cytoskeletal")
  }
  if (grepl("metabol|dehydrogenase|synthetase|hydrolase|transferase|oxidoreductase|enzyme", txt)) {
    return("Metabolism")
  }
  if (grepl("ribosome|translation|rna processing|splice|elongation", txt)) {
    return("RNA / translation")
  }

  return("Uncharacterized")
}

# ------------------------------------------------------------
# Load UV gene-level table
# ------------------------------------------------------------
cat("Loading UV full gene-level table...\n")
uv <- fread(uv_file)

if ("feature" %in% names(uv) && !("gene" %in% names(uv))) {
  setnames(uv, "feature", "gene")
}
if (!("gene" %in% names(uv))) stop("UV table has no gene column")
if (!("delta_pool" %in% names(uv))) stop("UV table missing delta_pool")
if (!("instability_score" %in% names(uv))) stop("UV table missing instability_score")

uv <- uv[!is.na(gene) & gene != ""]
uv[, gene := normalize_gene(gene)]
uv[, gene_base := base_gene(gene)]

# ------------------------------------------------------------
# Define formal masked genes
# ------------------------------------------------------------
cat("Defining masked genes using |delta_pool| <= Q25 and instability_score >= Q90...\n")
q25_pool <- quantile(abs(uv$delta_pool), 0.25, na.rm = TRUE)
q90_inst <- quantile(uv$instability_score, 0.90, na.rm = TRUE)

masked <- uv[
  abs(delta_pool) <= q25_pool &
    instability_score >= q90_inst
]

masked <- masked %>%
  as_tibble() %>%
  distinct(gene, .keep_all = TRUE)

write_tsv(masked, masked_out)
cat("Masked genes:", nrow(masked), "\n")

# ------------------------------------------------------------
# Load reference map
# ------------------------------------------------------------
cat("Loading NCBI reference map...\n")
ref <- read_csv(ref_file, show_col_types = FALSE) %>%
  mutate(
    gene       = normalize_gene(gene),
    locus_tag  = normalize_gene(locus_tag),
    gene_base  = base_gene(gene),
    locus_base = base_gene(locus_tag)
  )

# ------------------------------------------------------------
# Match masked genes to reference
# ------------------------------------------------------------
cat("Matching masked genes to NCBI reference...\n")

m_exact_gene <- masked %>%
  left_join(ref, by = c("gene" = "gene")) %>%
  filter(!is.na(protein_id)) %>%
  mutate(match_type = "exact_gene")

matched_1 <- unique(m_exact_gene$gene)

m_exact_locus <- masked %>%
  filter(!(gene %in% matched_1)) %>%
  left_join(ref, by = c("gene" = "locus_tag")) %>%
  filter(!is.na(protein_id)) %>%
  mutate(match_type = "exact_locus")

matched_2 <- unique(c(matched_1, m_exact_locus$gene))

m_base_gene <- masked %>%
  filter(!(gene %in% matched_2)) %>%
  left_join(ref, by = c("gene_base" = "gene_base")) %>%
  filter(!is.na(protein_id)) %>%
  mutate(match_type = "base_gene")

matched_3 <- unique(c(matched_2, m_base_gene$gene))

m_base_locus <- masked %>%
  filter(!(gene %in% matched_3)) %>%
  left_join(ref, by = c("gene_base" = "locus_base")) %>%
  filter(!is.na(protein_id)) %>%
  mutate(match_type = "base_locus")

matched <- bind_rows(
  m_exact_gene,
  m_exact_locus,
  m_base_gene,
  m_base_locus
) %>%
  distinct(gene, .keep_all = TRUE)

unmatched <- masked %>%
  filter(!(gene %in% matched$gene))

write_csv(matched, matched_out)
write_csv(unmatched, unmatched_out)

cat("Matched genes:", nrow(matched), "\n")
cat("Unmatched genes:", nrow(unmatched), "\n")

# ------------------------------------------------------------
# Load merged annotation from old Figure 5C
# ------------------------------------------------------------
cat("Loading merged annotation table...\n")
ann <- read_csv(annot_file, show_col_types = FALSE)

# inspect likely useful columns
candidate_text_cols <- c(
  "display_annot",
  "protein_desc_from_cds",
  "protein_desc_from_protein",
  "Description",
  "Preferred_name",
  "ipr_desc",
  "pfam_hits",
  "annot_text"
)

available_text_cols <- candidate_text_cols[candidate_text_cols %in% names(ann)]
if (length(available_text_cols) == 0) {
  stop("No usable annotation text columns found in merged annotation file.")
}

# harmonize gene column name if needed
if (!("gene" %in% names(ann))) {
  gene_candidates <- c("feature", "gene_id")
  found_gene <- gene_candidates[gene_candidates %in% names(ann)][1]
  if (is.na(found_gene)) stop("No gene-like column found in annotation file.")
  ann <- ann %>% rename(gene = !!sym(found_gene))
}

ann <- ann %>%
  mutate(
    gene = normalize_gene(gene),
    annot_text = pmap_chr(
      select(., any_of(available_text_cols)),
      ~ paste(clean_txt(c(...)), collapse = " ; ")
    )
  )

# ------------------------------------------------------------
# Merge annotation onto matched genes
# ------------------------------------------------------------
cat("Merging matched genes with annotation...\n")

annot_merged <- matched %>%
  left_join(ann, by = "gene") %>%
  mutate(
    annot_text = clean_txt(annot_text),
    category = vapply(annot_text, assign_category, character(1))
  )

if (nrow(unmatched) > 0) {
  unmatched_add <- unmatched %>%
    mutate(
      protein_id = NA_character_,
      match_type = "unmatched",
      annot_text = "",
      category = "Uncharacterized"
    )
  annot_merged <- bind_rows(annot_merged, unmatched_add)
}

write_csv(annot_merged, annot_merged_out)

# ------------------------------------------------------------
# Annotation coverage
# ------------------------------------------------------------
coverage <- annot_merged %>%
  mutate(has_annotation = annot_text != "") %>%
  summarise(
    total = n(),
    annotated = sum(has_annotation),
    unannotated = sum(!has_annotation),
    annotated_pct = round(100 * annotated / total, 1)
  )

write_tsv(coverage, coverage_out)

cat("Annotation coverage:\n")
print(coverage)

# ------------------------------------------------------------
# Count categories
# ------------------------------------------------------------
cat("Counting categories...\n")

counts <- annot_merged %>%
  count(category, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

write_tsv(counts, counts_out)

print(counts)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
cat("Making Figure 5E plot...\n")

p <- ggplot(counts, aes(x = reorder(category, n), y = n)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "Functional structure of instability-recovered masked genes",
    x = "Functional category",
    y = "Number of genes"
  ) +
  theme_classic(base_size = 12)

ggsave(plot_png, p, width = 7, height = 5, dpi = 300)
ggsave(plot_pdf, p, width = 7, height = 5)

cat("✅ Figure 5E outputs written to:\n")
cat("   ", work_dir, "\n")
