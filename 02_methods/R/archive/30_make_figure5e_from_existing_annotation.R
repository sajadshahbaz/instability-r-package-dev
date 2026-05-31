#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT not set")

# -----------------------------
# Paths
# -----------------------------
ann_dir <- file.path(
  nm2_root,
  "04_results", "old", "figure5c_ncbi_clean", "annotation"
)

eggnog_file   <- file.path(ann_dir, "figure5c_eggnog.emapper.annotations")
interpro_file <- file.path(ann_dir, "figure5c_interproscan.tsv")
pfam_file     <- file.path(ann_dir, "figure5c_pfam.domtblout")

map_file <- file.path(
  nm2_root,
  "04_results", "figure5e_global_functional_structure",
  "annotation_from_scratch",
  "figure5e_masked_gene_protein_table.tsv"
)

out_dir <- file.path(
  nm2_root,
  "04_results", "figure5e_global_functional_structure_final"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
clean_txt <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

collapse_unique <- function(x) {
  x <- clean_txt(x)
  x <- x[x != ""]
  if (!length(x)) return("")
  paste(unique(x), collapse = "; ")
}

read_eggnog_robust <- function(fp, proteins) {
  empty <- data.frame(
    protein_id = character(),
    eggnog_desc = character(),
    eggnog_pref = character(),
    stringsAsFactors = FALSE
  )

  if (!file.exists(fp)) return(empty)

  lines <- readLines(fp, warn = FALSE)
  lines_data <- lines[!startsWith(lines, "#")]

  if (!length(lines_data)) return(empty)

  first <- strsplit(lines_data[1], "\t")[[1]]

  # Case 1: header exists
  if (any(grepl("query", first, ignore.case = TRUE))) {
    x <- suppressWarnings(read_tsv(fp, comment = "#", show_col_types = FALSE))
    nms <- names(x)
    qcol <- nms[grepl("^query$|^#query$", nms, ignore.case = TRUE)][1]
    dcol <- nms[grepl("^description$", nms, ignore.case = TRUE)][1]
    pcol <- nms[grepl("preferred_name", nms, ignore.case = TRUE)][1]

    if (is.na(qcol)) return(empty)

    out <- x %>%
      transmute(
        protein_id = as.character(.data[[qcol]]),
        eggnog_desc = if (!is.na(dcol)) clean_txt(.data[[dcol]]) else "",
        eggnog_pref = if (!is.na(pcol)) clean_txt(.data[[pcol]]) else ""
      ) %>%
      filter(protein_id %in% proteins)

    if (nrow(out) == 0) return(empty)
    return(as.data.frame(out))
  }

  # Case 2: no header
  x <- suppressWarnings(read_tsv(fp, comment = "#", col_names = FALSE, show_col_types = FALSE))
  if (nrow(x) == 0 || ncol(x) < 1) return(empty)

  desc_col <- if (ncol(x) >= 8) 8 else NA_integer_
  pref_col <- if (ncol(x) >= 9) 9 else NA_integer_

  out <- tibble(
    protein_id = clean_txt(x[[1]]),
    eggnog_desc = if (!is.na(desc_col)) clean_txt(x[[desc_col]]) else "",
    eggnog_pref = if (!is.na(pref_col)) clean_txt(x[[pref_col]]) else ""
  ) %>%
    filter(protein_id %in% proteins)

  if (nrow(out) == 0) return(empty)
  as.data.frame(out)
}

read_interpro_robust <- function(fp, proteins) {
  empty <- data.frame(
    protein_id = character(),
    interpro_sig_desc = character(),
    interpro_desc = character(),
    stringsAsFactors = FALSE
  )

  if (!file.exists(fp)) return(empty)

  x <- suppressWarnings(read_tsv(fp, col_names = FALSE, show_col_types = FALSE))
  if (nrow(x) == 0) return(empty)

  out <- tibble(
    protein_id = clean_txt(x[[1]]),
    interpro_sig_desc = if (ncol(x) >= 6) clean_txt(x[[6]]) else "",
    interpro_desc = if (ncol(x) >= 13) clean_txt(x[[13]]) else ""
  ) %>%
    filter(protein_id %in% proteins) %>%
    group_by(protein_id) %>%
    summarise(
      interpro_sig_desc = collapse_unique(interpro_sig_desc),
      interpro_desc = collapse_unique(interpro_desc),
      .groups = "drop"
    )

  if (nrow(out) == 0) return(empty)
  as.data.frame(out)
}

read_pfam_robust <- function(fp, proteins) {
  empty <- data.frame(
    protein_id = character(),
    pfam_desc = character(),
    stringsAsFactors = FALSE
  )

  if (!file.exists(fp)) return(empty)

  lines <- readLines(fp, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  if (!length(lines)) return(empty)

  split_lines <- strsplit(lines, "\\s+")
  df <- lapply(split_lines, function(z) {
    if (length(z) < 4) return(NULL)
    data.frame(
      pfam_acc = z[1],
      protein_id = z[4],
      stringsAsFactors = FALSE
    )
  })
  df <- bind_rows(df)
  if (nrow(df) == 0) return(empty)

  out <- df %>%
    filter(protein_id %in% proteins) %>%
    group_by(protein_id) %>%
    summarise(
      pfam_desc = collapse_unique(pfam_acc),
      .groups = "drop"
    )

  if (nrow(out) == 0) return(empty)
  as.data.frame(out)
}

assign_category <- function(txt) {
  txt <- tolower(clean_txt(txt))
  if (txt == "") return("Uncharacterized")

  if (str_detect(txt, "gpcr|receptor|g-protein|kinase|phosphatase|signal")) {
    return("Signaling/Receptor")
  }
  if (str_detect(txt, "transcription|dna-binding|rna-binding|regulator|chromatin")) {
    return("Regulation")
  }
  if (str_detect(txt, "transport|channel|carrier|transporter|membrane")) {
    return("Transport/Membrane")
  }
  if (str_detect(txt, "protease|peptidase|proteolysis|ubiquitin|trypsin|serine protease")) {
    return("Proteolysis")
  }
  if (str_detect(txt, "stress|heat shock|oxidative|desiccation")) {
    return("Stress response")
  }
  if (str_detect(txt, "repair|replication|recombination")) {
    return("DNA repair / genome maintenance")
  }
  if (str_detect(txt, "cytoskeleton|actin|tubulin|structural|calycin")) {
    return("Structural / cytoskeletal")
  }
  if (str_detect(txt, "metabol|dehydrogenase|synthetase|hydrolase|transferase|oxidoreductase|enzyme")) {
    return("Metabolism")
  }
  if (str_detect(txt, "ribosome|translation|splice|elongation|rna processing")) {
    return("RNA / translation")
  }

  return("Other")
}

# -----------------------------
# Load masked proteins
# -----------------------------
mp <- fread(map_file)
proteins <- unique(mp$protein_id)

cat("Masked proteins:", length(proteins), "\n")
cat("Example protein IDs from masked set:\n")
print(head(proteins, 5))

# -----------------------------
# eggNOG diagnostics
# -----------------------------
if (file.exists(eggnog_file)) {
  raw_lines <- readLines(eggnog_file, warn = FALSE)
  raw_lines <- raw_lines[!startsWith(raw_lines, "#")]
  if (length(raw_lines)) {
    egg_ids <- sapply(strsplit(raw_lines[1:min(5, length(raw_lines))], "\t"), `[`, 1)
    cat("Example eggNOG IDs:\n")
    print(egg_ids)
  }
}

# -----------------------------
# Read annotation files
# -----------------------------
cat("Reading eggNOG...\n")
eggnog <- read_eggnog_robust(eggnog_file, proteins)
cat("eggNOG matched proteins:", n_distinct(eggnog$protein_id), "\n")

cat("Reading InterPro...\n")
interpro <- read_interpro_robust(interpro_file, proteins)
cat("InterPro matched proteins:", n_distinct(interpro$protein_id), "\n")

cat("Reading Pfam...\n")
pfam <- read_pfam_robust(pfam_file, proteins)
cat("Pfam matched proteins:", n_distinct(pfam$protein_id), "\n")

# -----------------------------
# Merge
# -----------------------------
ann <- mp %>%
  left_join(eggnog, by = "protein_id") %>%
  left_join(interpro, by = "protein_id") %>%
  left_join(pfam, by = "protein_id") %>%
  mutate(
    annotation_text = paste(
      clean_txt(eggnog_pref),
      clean_txt(eggnog_desc),
      clean_txt(interpro_sig_desc),
      clean_txt(interpro_desc),
      clean_txt(pfam_desc),
      sep = " ; "
    ),
    has_annotation = annotation_text != "",
    category = vapply(annotation_text, assign_category, character(1))
  )

write_tsv(ann, file.path(out_dir, "figure5e_annotation_merged.tsv"))

coverage <- ann %>%
  summarise(
    total = n(),
    annotated = sum(has_annotation),
    unannotated = sum(!has_annotation),
    annotated_pct = round(100 * annotated / total, 1)
  )

write_tsv(coverage, file.path(out_dir, "figure5e_annotation_coverage.tsv"))

cat("\nAnnotation coverage:\n")
print(coverage)

summary_table <- ann %>%
  count(category, sort = TRUE)

write_tsv(summary_table, file.path(out_dir, "figure5e_category_counts.tsv"))

cat("\nCategory counts:\n")
print(summary_table)

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(summary_table, aes(x = reorder(category, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Functional structure of instability-recovered masked genes",
    x = "Category",
    y = "Number of genes"
  )

ggsave(file.path(out_dir, "figure5e.png"), p, width = 6, height = 4)
ggsave(file.path(out_dir, "figure5e.pdf"), p, width = 6, height = 4)

cat("\n✅ Figure 5E generated in:\n", out_dir, "\n")
