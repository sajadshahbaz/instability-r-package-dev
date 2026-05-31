#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
uv_file <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")

ncbi_dir <- file.path(
  nm2_root,
  "01_data", "raw", "ncbi_dataset", "ncbi_dataset", "data", "GCA_001949185.1"
)
cds_fna <- file.path(ncbi_dir, "cds_from_genomic.fna")

ann_dir <- file.path(
  nm2_root,
  "04_results", "figure5e_global_functional_structure", "annotation_from_scratch"
)
map_file  <- file.path(ann_dir, "figure5e_masked_gene_protein_table.tsv")
pfam_file <- file.path(ann_dir, "figure5e_pfam.domtblout")
ipr_file  <- file.path(ann_dir, "figure5e_interproscan.tsv")

out_dir <- file.path(nm2_root, "04_results", "figure5e_final_refined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

masked_out         <- file.path(out_dir, "figure5e_masked_genes.tsv")
gene_map_out       <- file.path(out_dir, "figure5e_gene_protein_mapping.tsv")
annot_tbl_out      <- file.path(out_dir, "figure5e_annotation_table.tsv")
coverage_out       <- file.path(out_dir, "figure5e_annotation_coverage.tsv")
counts_out         <- file.path(out_dir, "figure5e_category_counts.tsv")
counts_refined_out <- file.path(out_dir, "figure5e_category_counts_refined.tsv")
plot_png           <- file.path(out_dir, "figure5e_refined.png")
plot_pdf           <- file.path(out_dir, "figure5e_refined.pdf")

for (fp in c(uv_file, cds_fna, map_file, pfam_file, ipr_file)) {
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

collapse_unique <- function(x) {
  x <- clean_txt(x)
  x <- x[x != "" & x != "-"]
  if (!length(x)) return("")
  paste(unique(x), collapse = "; ")
}

extract_tag <- function(header, tag) {
  m <- str_match(header, paste0("\\[", tag, "=([^\\]]+)\\]"))
  m[, 2]
}

assign_category <- function(txt) {
  txt <- tolower(clean_txt(txt))
  txt <- gsub("\\s+", " ", txt)

  if (txt == "") return("Uncharacterized")

  if (grepl("receptor|gpcr|kinase|phosphatase|g-protein|signal|pdz|lrr|leucine rich repeat", txt)) {
    return("Signaling/Receptor")
  }
  if (grepl("protease|peptidase|trypsin|serine protease|cathepsin|metalloprotease", txt)) {
    return("Proteolysis")
  }
  if (grepl("transporter|channel|carrier|membrane|slc|sulp|stas", txt)) {
    return("Transport/Membrane")
  }
  if (grepl("transcription|dna-binding|rna-binding|regulator|chromatin|zinc finger", txt)) {
    return("Regulation/Nucleic-acid associated")
  }
  if (grepl("dehydrogenase|synthetase|hydrolase|transferase|oxidoreductase|enzyme|metabol", txt)) {
    return("Metabolism/Enzyme")
  }
  if (grepl("actin|tubulin|cytoskeleton|structural|motor protein|coiled-coil", txt)) {
    return("Structural/Other")
  }
  if (grepl("ribosome|translation|splice|elongation|rna processing", txt)) {
    return("RNA / translation")
  }

  return("Other annotated")
}

refine_other <- function(txt) {
  txt <- tolower(clean_txt(txt))
  txt <- gsub("\\s+", " ", txt)

  if (txt == "") return("Uncharacterized")
  if (grepl("duf", txt)) return("DUF / unknown conserved")
  if (grepl("repeat|lrr|leucine rich repeat|wd40|ankyrin|tpr|armadillo|helical repeat", txt)) {
    return("Repeat / scaffold proteins")
  }
  if (grepl("zinc finger|dna-binding|rna-binding|binding protein", txt)) {
    return("Binding proteins")
  }
  if (grepl("coiled-coil|structural", txt)) {
    return("Coiled-coil / structural")
  }
  if (grepl("domain", txt)) {
    return("Generic domain proteins")
  }

  return("Other complex proteins")
}

# ------------------------------------------------------------
# 1) Formal masked set
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
# 2) Parse CDS headers to recover gene/locus/protein_id
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
# 3) Map masked genes to proteins
# ------------------------------------------------------------
cat("Mapping masked genes to protein IDs...\n")

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

mapped <- bind_rows(matched1, m2) %>%
  distinct(gene, .keep_all = TRUE) %>%
  select(gene, gene_base, protein_id, protein_desc, match_type)

write_tsv(mapped, gene_map_out)
cat("Mapped proteins:", nrow(mapped), "\n")

# ------------------------------------------------------------
# 4) Parse fresh Pfam domtblout
# ------------------------------------------------------------
cat("Parsing Pfam output...\n")
pfam_lines <- readLines(pfam_file, warn = FALSE)
pfam_lines <- pfam_lines[!startsWith(pfam_lines, "#")]

pfam_df <- lapply(pfam_lines, function(ln) {
  z <- strsplit(trimws(ln), "\\s+")[[1]]
  if (length(z) < 23) return(NULL)
  data.frame(
    pfam_name  = z[1],
    pfam_acc   = z[2],
    protein_id = z[4],
    pfam_desc  = paste(z[23:length(z)], collapse = " "),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

pfam_sum <- pfam_df %>%
  group_by(protein_id) %>%
  summarise(
    pfam_name = collapse_unique(pfam_name),
    pfam_acc  = collapse_unique(pfam_acc),
    pfam_desc = collapse_unique(pfam_desc),
    .groups = "drop"
  )

cat("Pfam-annotated proteins:", nrow(pfam_sum), "\n")

# ------------------------------------------------------------
# 5) Parse fresh InterProScan TSV
# ------------------------------------------------------------
cat("Parsing InterProScan output...\n")
ipr <- read_tsv(ipr_file, col_names = FALSE, show_col_types = FALSE)

ipr_sum <- ipr %>%
  transmute(
    protein_id = clean_txt(X1),
    analysis   = clean_txt(X4),
    sig_acc    = clean_txt(X5),
    sig_desc   = clean_txt(X6),
    ipr_acc    = if ("X12" %in% names(.)) clean_txt(X12) else "",
    ipr_desc   = if ("X13" %in% names(.)) clean_txt(X13) else ""
  ) %>%
  group_by(protein_id) %>%
  summarise(
    analysis = collapse_unique(analysis),
    sig_acc  = collapse_unique(sig_acc),
    sig_desc = collapse_unique(sig_desc),
    ipr_acc  = collapse_unique(ipr_acc),
    ipr_desc = collapse_unique(ipr_desc),
    .groups = "drop"
  )

cat("InterPro-annotated proteins:", nrow(ipr_sum), "\n")

# ------------------------------------------------------------
# 6) Merge and classify
# ------------------------------------------------------------
cat("Merging annotations and assigning categories...\n")
annot_tbl <- mapped %>%
  left_join(pfam_sum, by = "protein_id") %>%
  left_join(ipr_sum, by = "protein_id") %>%
  rowwise() %>%
  mutate(
    annotation_text = collapse_unique(c(
      protein_desc,
      pfam_desc,
      sig_desc,
      ipr_desc,
      pfam_name,
      analysis
    )),
    has_annotation = annotation_text != "",
    category = ifelse(has_annotation, assign_category(annotation_text), "Uncharacterized")
  ) %>%
  ungroup() %>%
  mutate(
    category_refined = case_when(
      category != "Other annotated" ~ category,
      TRUE ~ vapply(annotation_text, refine_other, character(1))
    )
  )

write_tsv(annot_tbl, annot_tbl_out)

# ------------------------------------------------------------
# 7) Coverage
# ------------------------------------------------------------
coverage <- annot_tbl %>%
  summarise(
    total = n(),
    annotated = sum(has_annotation),
    unannotated = sum(!has_annotation),
    annotated_pct = round(100 * annotated / total, 1)
  )

write_tsv(coverage, coverage_out)

cat("\nAnnotation coverage:\n")
print(coverage)

# ------------------------------------------------------------
# 8) Counts
# ------------------------------------------------------------
counts <- annot_tbl %>%
  count(category, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

counts_refined <- annot_tbl %>%
  count(category_refined, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

write_tsv(counts, counts_out)
write_tsv(counts_refined, counts_refined_out)

cat("\nBase category counts:\n")
print(counts)

cat("\nRefined category counts:\n")
print(counts_refined)

# ------------------------------------------------------------
# 9) Plot refined
# ------------------------------------------------------------
p <- ggplot(counts_refined, aes(x = reorder(category_refined, n), y = n)) +
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

cat("\n✅ Final Figure 5E outputs written to:\n", out_dir, "\n")
