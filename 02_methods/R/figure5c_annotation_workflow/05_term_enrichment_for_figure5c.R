#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(forcats)
})

# ============================================================
# Figure 5C — Term/domain enrichment for Instability vs Pooled
# ============================================================

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT is not set")
work_dir <- file.path(nm2_root, "04_results", "figure5c_ncbi_clean")
fig_dir  <- file.path(work_dir, "figure5c_final")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

infile <- file.path(fig_dir, "figure5c_fresh_annotation_merged.csv")
if (!file.exists(infile)) stop("Missing merged annotation file: ", infile)

df <- read_csv(infile, show_col_types = FALSE) %>%
  filter(category %in% c("Instability-only", "Pooled-only")) %>%
  distinct(gene, category, .keep_all = TRUE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
split_terms <- function(x, sep = ";") {
  x <- as.character(x)
  x[is.na(x)] <- ""
  strsplit(x, sep, fixed = TRUE)
}

make_long_terms <- function(df, col, label_col = NULL, source_name = NULL) {
  vals <- split_terms(df[[col]], ";")
  out <- purrr::map2_dfr(vals, seq_along(vals), function(v, i) {
    v <- trimws(v)
    v <- v[v != "" & !is.na(v)]
    if (!length(v)) return(NULL)
    tibble(
      gene = df$gene[i],
      category = df$category[i],
      term = v
    )
  })
  if (nrow(out) == 0) return(tibble(gene=character(), category=character(), term=character(), source=character()))
  out %>%
    mutate(source = ifelse(is.null(source_name), col, source_name)) %>%
    distinct()
}

fisher_enrich <- function(long_df, source_name, min_total = 3) {
  if (nrow(long_df) == 0) return(tibble())

  genes_inst <- unique(long_df$gene[long_df$category == "Instability-only"])
  genes_pool <- unique(long_df$gene[long_df$category == "Pooled-only"])

  n_inst <- length(unique(df$gene[df$category == "Instability-only"]))
  n_pool <- length(unique(df$gene[df$category == "Pooled-only"]))

  term_tbl <- long_df %>%
    group_by(term) %>%
    summarise(
      inst_hits = n_distinct(gene[category == "Instability-only"]),
      pool_hits = n_distinct(gene[category == "Pooled-only"]),
      total_hits = n_distinct(gene),
      .groups = "drop"
    ) %>%
    filter(total_hits >= min_total)

  if (nrow(term_tbl) == 0) return(tibble())

  res <- term_tbl %>%
    rowwise() %>%
    mutate(
      inst_no = n_inst - inst_hits,
      pool_no = n_pool - pool_hits,
      fisher_p = fisher.test(matrix(c(inst_hits, inst_no, pool_hits, pool_no), nrow = 2))$p.value,
      odds_ratio = unname(fisher.test(matrix(c(inst_hits, inst_no, pool_hits, pool_no), nrow = 2))$estimate)
    ) %>%
    ungroup() %>%
    mutate(
      padj = p.adjust(fisher_p, method = "BH"),
      log2_or = log2(odds_ratio),
      enriched_in = case_when(
        is.infinite(log2_or) & inst_hits > pool_hits ~ "Instability-only",
        is.infinite(log2_or) & pool_hits > inst_hits ~ "Pooled-only",
        log2_or > 0 ~ "Instability-only",
        log2_or < 0 ~ "Pooled-only",
        TRUE ~ "Neither"
      ),
      source = source_name
    ) %>%
    arrange(padj, desc(abs(log2_or)))

  res
}

# ------------------------------------------------------------
# Build term tables
# ------------------------------------------------------------
pfam_long <- make_long_terms(df, "pfam_hits", source_name = "Pfam")
ipracc_long <- make_long_terms(df, "ipr_acc", source_name = "InterPro accession")
iprdesc_long <- make_long_terms(df, "ipr_desc", source_name = "InterPro description")
iprgo_long <- make_long_terms(df, "ipr_go", source_name = "InterPro GO")

# ------------------------------------------------------------
# Run enrichment
# ------------------------------------------------------------
pfam_res   <- fisher_enrich(pfam_long, "Pfam", min_total = 3)
ipracc_res <- fisher_enrich(ipracc_long, "InterPro accession", min_total = 3)
iprdesc_res <- fisher_enrich(iprdesc_long, "InterPro description", min_total = 3)
iprgo_res  <- fisher_enrich(iprgo_long, "InterPro GO", min_total = 3)

write_csv(pfam_res, file.path(fig_dir, "figure5c_pfam_enrichment.csv"))
write_csv(ipracc_res, file.path(fig_dir, "figure5c_interpro_accession_enrichment.csv"))
write_csv(iprdesc_res, file.path(fig_dir, "figure5c_interpro_description_enrichment.csv"))
write_csv(iprgo_res, file.path(fig_dir, "figure5c_interpro_go_enrichment.csv"))

# ------------------------------------------------------------
# Combine for plotting
# Prefer description-level because it is more readable
# ------------------------------------------------------------
plot_candidates <- bind_rows(
  pfam_res %>% mutate(display_term = term),
  iprdesc_res %>% mutate(display_term = term)
) %>%
  filter(!is.na(log2_or), !is.na(padj)) %>%
  filter(term != "") %>%
  filter(total_hits >= 3)

# choose top informative terms
plot_df <- plot_candidates %>%
  filter(padj <= 0.2 | abs(log2_or) >= 1) %>%
  arrange(padj, desc(abs(log2_or))) %>%
  slice_head(n = 20) %>%
  mutate(
    display_term = str_trunc(display_term, 70),
    display_term = fct_reorder(display_term, log2_or)
  )

write_csv(plot_df, file.path(fig_dir, "figure5c_plot_terms.csv"))

if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df, aes(x = log2_or, y = display_term, color = enriched_in, size = -log10(padj))) +
    geom_point(alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Term/domain enrichment in UV method-specific genes",
      x = "log2(odds ratio), Instability-only vs Pooled-only",
      y = NULL,
      color = NULL,
      size = "-log10(FDR)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(fig_dir, "Figure5C_term_enrichment.pdf"), p, width = 9, height = 6)
  ggsave(file.path(fig_dir, "Figure5C_term_enrichment.png"), p, width = 9, height = 6, dpi = 300)
}

# ------------------------------------------------------------
# Representative annotated genes
# ------------------------------------------------------------
repr <- df %>%
  mutate(
    best_annotation = coalesce(ipr_desc, pfam_hits, protein_desc_from_protein, protein_desc_from_cds)
  ) %>%
  filter(!is.na(best_annotation), best_annotation != "") %>%
  group_by(category) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  select(gene, category, protein_id, protein_desc_from_protein, ipr_desc, pfam_hits)

write_csv(repr, file.path(fig_dir, "figure5c_representative_annotated_genes.csv"))

cat("✅ Wrote enrichment tables and representative-gene table to:\n", fig_dir, "\n")
cat("Pfam terms tested:", nrow(pfam_res), "\n")
cat("InterPro description terms tested:", nrow(iprdesc_res), "\n")
cat("Plot terms retained:", nrow(plot_df), "\n")
