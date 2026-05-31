#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# ============================================================
# Figure 5C — Merge fresh annotation and plot
# Robust version: eggNOG may be empty
# ============================================================

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT is not set")
work_dir <- file.path(nm2_root, "04_results", "figure5c_ncbi_clean")
ann_dir <- file.path(work_dir, "annotation")
out_dir <- file.path(work_dir, "figure5c_final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

match_file  <- file.path(work_dir, "figure5c_genes_matched_to_ncbi.csv")
eggnog_file <- file.path(ann_dir, "figure5c_eggnog.emapper.annotations")
ipr_file    <- file.path(ann_dir, "figure5c_interproscan.tsv")
pfam_file   <- file.path(ann_dir, "figure5c_pfam.domtblout")

if (!file.exists(match_file)) stop("Missing file: ", match_file)

# -----------------------------
# Load matched genes
# -----------------------------
base_df <- read_csv(match_file, show_col_types = FALSE) %>%
  distinct(gene, category, protein_id, protein_desc_from_cds, protein_desc_from_protein)

# -----------------------------
# eggNOG parser
# -----------------------------
read_eggnog <- function(fp) {
  empty <- tibble(
    protein_id = character(),
    Preferred_name = character(),
    Description = character(),
    GOs = character(),
    KEGG_ko = character(),
    KEGG_Pathway = character(),
    PFAMs_eggnog = character()
  )

  if (!file.exists(fp)) return(empty)

  lines <- read_lines(fp)
  lines <- lines[!startsWith(lines, "#")]
  if (!length(lines)) return(empty)

  hdr <- strsplit(lines[1], "\t")[[1]]
  has_header <- any(grepl("query", hdr, ignore.case = TRUE))

  df <- tryCatch({
    if (has_header) {
      read_tsv(fp, comment = "#", show_col_types = FALSE)
    } else {
      col_names <- c(
        "query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl",
        "COG_cat","Description","Preferred_name","GOs","EC","KEGG_ko",
        "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
        "BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"
      )
      suppressWarnings(read_tsv(fp, comment = "#", col_names = col_names, show_col_types = FALSE))
    }
  }, error = function(e) NULL)

  if (is.null(df) || nrow(df) == 0) return(empty)

  out <- df %>%
    transmute(
      protein_id = query,
      Preferred_name = coalesce(Preferred_name, ""),
      Description = coalesce(Description, ""),
      GOs = coalesce(GOs, ""),
      KEGG_ko = coalesce(KEGG_ko, ""),
      KEGG_Pathway = coalesce(KEGG_Pathway, ""),
      PFAMs_eggnog = coalesce(PFAMs, "")
    ) %>%
    filter(!is.na(protein_id), protein_id != "")

  if (nrow(out) == 0) return(empty)
  out
}

# -----------------------------
# InterPro parser
# -----------------------------
read_interpro <- function(fp) {
  empty <- tibble(
    protein_id = character(),
    ipr_acc = character(),
    ipr_desc = character(),
    ipr_go = character(),
    ipr_path = character()
  )

  if (!file.exists(fp)) return(empty)

  x <- tryCatch(read_tsv(fp, col_names = FALSE, show_col_types = FALSE),
                error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(empty)

  ncols <- ncol(x)
  colnames(x) <- paste0("V", seq_len(ncols))

  getcol <- function(df, nm) {
    if (nm %in% colnames(df)) as.character(df[[nm]]) else rep(NA_character_, nrow(df))
  }

  out <- tibble(
    protein_id = getcol(x, "V1"),
    ipr_acc    = getcol(x, "V12"),
    ipr_desc   = getcol(x, "V13"),
    ipr_go     = getcol(x, "V14"),
    ipr_path   = getcol(x, "V15")
  ) %>%
    filter(!is.na(protein_id), protein_id != "") %>%
    group_by(protein_id) %>%
    summarise(
      ipr_acc  = paste(unique(na.omit(ipr_acc)), collapse = ";"),
      ipr_desc = paste(unique(na.omit(ipr_desc)), collapse = ";"),
      ipr_go   = paste(unique(na.omit(ipr_go)), collapse = ";"),
      ipr_path = paste(unique(na.omit(ipr_path)), collapse = ";"),
      .groups = "drop"
    )

  if (nrow(out) == 0) return(empty)
  out
}

# -----------------------------
# Pfam parser
# -----------------------------
read_pfam <- function(fp) {
  empty <- tibble(
    protein_id = character(),
    pfam_hits = character()
  )

  if (!file.exists(fp)) return(empty)

  lines <- read_lines(fp)
  lines <- lines[!startsWith(lines, "#")]
  if (!length(lines)) return(empty)

  parts <- strsplit(lines, "\\s+")
  df <- purrr::map_dfr(parts, function(x) {
    if (length(x) < 4) return(NULL)
    tibble(
      pfam_id = x[1],
      protein_id = x[4]
    )
  })

  if (nrow(df) == 0) return(empty)

  out <- df %>%
    filter(!is.na(protein_id), protein_id != "") %>%
    group_by(protein_id) %>%
    summarise(
      pfam_hits = paste(unique(pfam_id), collapse = ";"),
      .groups = "drop"
    )

  if (nrow(out) == 0) return(empty)
  out
}

# -----------------------------
# Read annotation sources
# -----------------------------
eg <- read_eggnog(eggnog_file)
ip <- read_interpro(ipr_file)
pf <- read_pfam(pfam_file)

cat("eggNOG rows:", nrow(eg), " cols:", paste(colnames(eg), collapse = ", "), "\n")
cat("InterPro rows:", nrow(ip), " cols:", paste(colnames(ip), collapse = ", "), "\n")
cat("Pfam rows:", nrow(pf), " cols:", paste(colnames(pf), collapse = ", "), "\n")

# -----------------------------
# Merge all annotations
# -----------------------------
merged <- base_df %>%
  left_join(eg, by = "protein_id") %>%
  left_join(ip, by = "protein_id") %>%
  left_join(pf, by = "protein_id")

write_csv(merged, file.path(out_dir, "figure5c_fresh_annotation_merged.csv"))

# -----------------------------
# Annotation coverage summary
# -----------------------------
coverage <- merged %>%
  mutate(
    has_ncbi_desc = (!is.na(protein_desc_from_cds) & protein_desc_from_cds != "") |
                    (!is.na(protein_desc_from_protein) & protein_desc_from_protein != ""),
    has_eggnog = (!is.na(Description) & Description != "") |
                 (!is.na(Preferred_name) & Preferred_name != "") |
                 (!is.na(GOs) & GOs != "") |
                 (!is.na(KEGG_Pathway) & KEGG_Pathway != ""),
    has_interpro = (!is.na(ipr_acc) & ipr_acc != "") |
                   (!is.na(ipr_desc) & ipr_desc != ""),
    has_pfam = (!is.na(pfam_hits) & pfam_hits != "")
  ) %>%
  group_by(category) %>%
  summarise(
    total = n(),
    ncbi_desc = sum(has_ncbi_desc),
    eggnog = sum(has_eggnog),
    interpro = sum(has_interpro),
    pfam = sum(has_pfam),
    any_annotation = sum(has_ncbi_desc | has_eggnog | has_interpro | has_pfam),
    any_annotation_pct = round(100 * any_annotation / total, 1),
    .groups = "drop"
  )

write_csv(coverage, file.path(out_dir, "figure5c_annotation_coverage.csv"))

# -----------------------------
# Function classes
# -----------------------------
classify_one <- function(desc, pref, go, kegg, pfam1, pfam2, ipr, fallback1, fallback2) {
  vals <- c(desc, pref, go, kegg, pfam1, pfam2, ipr, fallback1, fallback2)
  vals <- vals[!is.na(vals) & vals != ""]
  txt <- tolower(paste(vals, collapse = " "))

  if (txt == "") {
    return("Unknown / uncharacterized")
  } else if (str_detect(txt, "dna repair|photolyase|rad51|rad52|genome integrity|nucleotide excision|base excision|mismatch repair")) {
    return("DNA repair / genome integrity")
  } else if (str_detect(txt, "transcription|chromatin|histone|rna|splicing|nucleosome|transcription factor")) {
    return("Transcription / chromatin / RNA")
  } else if (str_detect(txt, "kinase|signal|receptor|gpcr|calcium|mapk|camp|g-protein")) {
    return("Signaling / regulation")
  } else if (str_detect(txt, "transporter|channel|membrane|vesicle|abc transporter|solute carrier")) {
    return("Transport / membrane")
  } else if (str_detect(txt, "chaperone|heat shock|proteasome|ubiquitin|autophagy|protein folding")) {
    return("Proteostasis / stress handling")
  } else if (str_detect(txt, "metabolism|oxidoreductase|redox|mitochond|glycolysis|biosynthesis|catabolism|dehydrogenase|synthetase|hydrolase")) {
    return("Metabolism / redox")
  } else if (str_detect(txt, "cytoskeleton|actin|tubulin|motor protein|myosin|dynein|kinesin")) {
    return("Structure / cytoskeleton")
  } else if (str_detect(txt, "unknown|uncharacterized|hypothetical")) {
    return("Unknown / uncharacterized")
  } else {
    return("Other annotated")
  }
}

annotated <- merged %>%
  mutate(
    function_class = purrr::pmap_chr(
      list(
        Description, Preferred_name, GOs, KEGG_Pathway,
        PFAMs_eggnog, pfam_hits, ipr_desc,
        protein_desc_from_cds, protein_desc_from_protein
      ),
      classify_one
    )
  )

write_csv(annotated, file.path(out_dir, "figure5c_annotated_with_classes.csv"))

# -----------------------------
# Summary and plot
# -----------------------------
sum_df <- annotated %>%
  group_by(category, function_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(category) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write_csv(sum_df, file.path(out_dir, "figure5c_summary.csv"))

keep_classes <- sum_df %>%
  group_by(function_class) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  filter(total_n >= 2) %>%
  arrange(total_n) %>%
  pull(function_class)

plot_df <- sum_df %>%
  filter(function_class %in% keep_classes) %>%
  mutate(
    category = factor(category, levels = c("Instability-only", "Pooled-only")),
    function_class = fct_reorder(function_class, n, .fun = sum)
  )

p <- ggplot(plot_df, aes(x = prop, y = function_class, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_x_continuous(labels = percent_format()) +
  labs(
    title = "Functional composition of method-specific genes (UV)",
    x = "Proportion within group",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "Figure5C_fresh.pdf"), p, width = 9, height = 5.5)
ggsave(file.path(out_dir, "Figure5C_fresh.png"), p, width = 9, height = 5.5, dpi = 300)

cat("✅ Wrote merged annotation, class table, coverage, summary, and Figure 5C to:\n", out_dir, "\n")
print(coverage)
