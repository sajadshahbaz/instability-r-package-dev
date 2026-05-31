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
base_dir <- file.path(
  nm2_root,
  "04_results", "figure5e_global_functional_structure",
  "annotation_from_scratch"
)

map_file  <- file.path(base_dir, "figure5e_masked_gene_protein_table.tsv")
pfam_file <- file.path(base_dir, "figure5e_pfam.domtblout")
ipr_file  <- file.path(base_dir, "figure5e_interproscan.tsv")

out_dir <- file.path(
  nm2_root,
  "04_results", "figure5e_fresh_annotation"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

annot_out    <- file.path(out_dir, "figure5e_fresh_annotation_table.tsv")
counts_out   <- file.path(out_dir, "figure5e_category_counts.tsv")
coverage_out <- file.path(out_dir, "figure5e_annotation_coverage.tsv")
plot_png     <- file.path(out_dir, "figure5e_fresh_annotation.png")
plot_pdf     <- file.path(out_dir, "figure5e_fresh_annotation.pdf")

for (fp in c(map_file, pfam_file, ipr_file)) {
  if (!file.exists(fp)) stop("Missing file: ", fp)
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
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

assign_category <- function(txt) {
  txt <- tolower(clean_txt(txt))
  txt <- gsub("\\s+", " ", txt)

  if (txt == "") return("Uncharacterized")

  if (grepl("leucine rich repeat|lrr|pdz|stas domain|kinase|phosphatase|g-protein|signal|receptor", txt)) {
    return("Signaling/Receptor")
  }
  if (grepl("protease|peptidase|trypsin|serine protease|cathepsin|metalloprotease", txt)) {
    return("Proteolysis")
  }
  if (grepl("transporter|channel|carrier|membrane|slc|sulp|stas", txt)) {
    return("Transport/Membrane")
  }
  if (grepl("transcription|dna-binding|rna-binding|regulator|chromatin", txt)) {
    return("Regulation/Nucleic-acid associated")
  }
  if (grepl("dehydrogenase|synthetase|hydrolase|transferase|oxidoreductase|enzyme|metabol", txt)) {
    return("Metabolism/Enzyme")
  }
  if (grepl("actin|tubulin|cytoskeleton|structural|motor protein", txt)) {
    return("Structural/Other")
  }
  if (grepl("ribosome|translation|splice|elongation|rna processing", txt)) {
    return("RNA / translation")
  }

  return("Other annotated")
}

# ------------------------------------------------------------
# Load gene-protein map
# ------------------------------------------------------------
mp <- read_tsv(map_file, show_col_types = FALSE) %>%
  distinct(gene, protein_id, .keep_all = TRUE)

cat("Mapped proteins:", n_distinct(mp$protein_id), "\n")

# ------------------------------------------------------------
# Parse Pfam domtblout
# ------------------------------------------------------------
pfam_lines <- readLines(pfam_file, warn = FALSE)
pfam_lines <- pfam_lines[!startsWith(pfam_lines, "#")]

pfam_df <- lapply(pfam_lines, function(ln) {
  z <- strsplit(trimws(ln), "\\s+")[[1]]
  if (length(z) < 23) return(NULL)
  data.frame(
    pfam_name = z[1],
    pfam_acc  = z[2],
    protein_id = z[4],
    pfam_desc = paste(z[23:length(z)], collapse = " "),
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
# Parse InterProScan TSV
# Columns:
# 1 protein_id
# 4 analysis
# 5 signature accession
# 6 signature desc
# 12 IPR accession
# 13 IPR desc
# ------------------------------------------------------------
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
# Merge
# ------------------------------------------------------------
annot_tbl <- mp %>%
  left_join(pfam_sum, by = "protein_id") %>%
  left_join(ipr_sum, by = "protein_id") %>%
  rowwise() %>%
  mutate(
    annotation_text = collapse_unique(c(
      pfam_desc,
      sig_desc,
      ipr_desc,
      pfam_name,
      analysis
    )),
    has_annotation = annotation_text != "",
    category = ifelse(has_annotation, assign_category(annotation_text), "Uncharacterized")
  ) %>%
  ungroup()

write_tsv(annot_tbl, annot_out)

# ------------------------------------------------------------
# Coverage
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
# Category counts
# ------------------------------------------------------------
counts <- annot_tbl %>%
  count(category, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

write_tsv(counts, counts_out)

cat("\nCategory counts:\n")
print(counts)

# ------------------------------------------------------------
# Plot
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
