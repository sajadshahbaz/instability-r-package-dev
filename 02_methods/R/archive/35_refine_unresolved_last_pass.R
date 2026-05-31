#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

in_dir <- file.path(nm2_root, "04_results", "figure5e_final_refined")
annot_file <- file.path(in_dir, "figure5e_annotation_table.tsv")

out_dir <- file.path(nm2_root, "04_results", "figure5e_final_refined_v3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

annot_out  <- file.path(out_dir, "figure5e_annotation_table_v3.tsv")
counts_out <- file.path(out_dir, "figure5e_category_counts_final_v3.tsv")
plot_png   <- file.path(out_dir, "figure5e_final_v3.png")
plot_pdf   <- file.path(out_dir, "figure5e_final_v3.pdf")
unresolved_out <- file.path(out_dir, "figure5e_unresolved_breakdown.tsv")

if (!file.exists(annot_file)) stop("Missing file: ", annot_file)

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

refine_unresolved <- function(txt) {
  txt <- tolower(clean_txt(txt))
  txt <- gsub("\\s+", " ", txt)

  if (txt == "") return("Unresolved: no detectable feature")

  # membrane / transporter / TM-like
  if (grepl("transmembrane|membrane|channel|transporter|permease|slc|sulp|stas", txt)) {
    return("Unresolved: membrane-associated")
  }

  # repeats / scaffold / domain-rich
  if (grepl("repeat|lrr|leucine rich repeat|ankyrin|wd40|armadillo|tpr|pdz|helical repeat|scaffold", txt)) {
    return("Unresolved: repeat/scaffold-rich")
  }

  # catalytic but broad / weakly specified
  if (grepl("hydrolase|transferase|oxidoreductase|synthetase|catalytic|enzyme|esterase|peptidase|nuclease", txt)) {
    return("Unresolved: catalytic-like")
  }

  # binding / regulatory-like
  if (grepl("binding|dna-binding|rna-binding|zinc finger|regulator|chromatin|adapter", txt)) {
    return("Unresolved: binding/regulatory-like")
  }

  # structural / coiled / low complexity / disordered
  if (grepl("coiled|structural|disorder|low complexity|mobidb|compositionally biased", txt)) {
    return("Unresolved: structural/disordered")
  }

  # conserved but generic
  if (grepl("conserved|superfamily|domain|family|motif|fold|putative|predicted", txt)) {
    return("Unresolved: conserved generic domain")
  }

  return("Unresolved: complex unknown")
}

ann <- read_tsv(annot_file, show_col_types = FALSE)

required_cols <- c("gene", "protein_id", "annotation_text", "category_refined")
miss <- setdiff(required_cols, names(ann))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

ann3 <- ann %>%
  mutate(
    category_final = case_when(
      category_refined != "Other complex proteins" ~ category_refined,
      TRUE ~ vapply(annotation_text, refine_unresolved, character(1))
    )
  )

write_tsv(ann3, annot_out)

unresolved_tbl <- ann3 %>%
  filter(str_detect(category_final, "^Unresolved:")) %>%
  count(category_final, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

write_tsv(unresolved_tbl, unresolved_out)

counts_final <- ann3 %>%
  count(category_final, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

write_tsv(counts_final, counts_out)

cat("\nUnresolved breakdown:\n")
print(unresolved_tbl)

cat("\nFinal category counts:\n")
print(counts_final)

p <- ggplot(counts_final, aes(x = reorder(category_final, n), y = n)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Functional structure of instability-recovered masked genes",
    x = "Broad functional class",
    y = "Number of genes"
  )

ggsave(plot_png, p, width = 7.2, height = 5.4)
ggsave(plot_pdf, p, width = 7.2, height = 5.4)

cat("\n✅ Final v3 outputs written to:\n", out_dir, "\n")
