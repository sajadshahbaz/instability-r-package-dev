#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
in_dir <- file.path(nm2_root, "04_results", "figure5e_final_refined")
annot_file <- file.path(in_dir, "figure5e_annotation_table.tsv")

out_dir <- file.path(nm2_root, "04_results", "figure5e_final_refined_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

annot_out  <- file.path(out_dir, "figure5e_annotation_table_v2.tsv")
counts_out <- file.path(out_dir, "figure5e_category_counts_final.tsv")
plot_png   <- file.path(out_dir, "figure5e_final_v2.png")
plot_pdf   <- file.path(out_dir, "figure5e_final_v2.pdf")

if (!file.exists(annot_file)) stop("Missing file: ", annot_file)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clean_txt <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

refine_other_complex <- function(txt) {
  txt <- tolower(clean_txt(txt))
  txt <- gsub("\\s+", " ", txt)

  if (txt == "") return("Complex proteins, unresolved class")

  if (grepl("transmembrane|membrane|channel|transporter|slc|sulp|stas", txt)) {
    return("Membrane-associated complex proteins")
  }

  if (grepl("hydrolase|transferase|oxidoreductase|synthetase|catalytic|enzyme|metallo", txt)) {
    return("Catalytic-like complex proteins")
  }

  if (grepl("repeat|motif|domain|superfamily|family|helical|fold", txt)) {
    return("Complex domain proteins")
  }

  if (grepl("conserved|putative|predicted", txt)) {
    return("Conserved complex proteins")
  }

  return("Complex proteins, unresolved class")
}

# ------------------------------------------------------------
# Load annotation table
# ------------------------------------------------------------
ann <- read_tsv(annot_file, show_col_types = FALSE)

required_cols <- c("gene", "protein_id", "annotation_text", "category_refined")
miss <- setdiff(required_cols, names(ann))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

# ------------------------------------------------------------
# Final refinement
# ------------------------------------------------------------
ann2 <- ann %>%
  mutate(
    category_final = case_when(
      category_refined != "Other complex proteins" ~ category_refined,
      TRUE ~ vapply(annotation_text, refine_other_complex, character(1))
    )
  )

write_tsv(ann2, annot_out)

counts_final <- ann2 %>%
  count(category_final, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

write_tsv(counts_final, counts_out)

cat("\nFinal category counts:\n")
print(counts_final)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(counts_final, aes(x = reorder(category_final, n), y = n)) +
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

cat("\n✅ Final refined outputs written to:\n", out_dir, "\n")
