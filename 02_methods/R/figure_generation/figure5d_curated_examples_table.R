#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(gridExtra)
  library(grid)
})

# ============================================================
# Figure 5D — Curated biological examples (TABLE, FINAL REVISED)
# End-to-end from raw Figure 5B sets + UV metrics + annotation
# ============================================================

nm2_root <- Sys.getenv("NM2_ROOT")
if (nm2_root == "") stop("NM2_ROOT is not set")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
fig5b_dir <- file.path(nm2_root, "04_results", "figure5b_topN_overlap")
uv_file   <- file.path(nm2_root, "04_results", "tables", "uv_pool_vs_strat_all.tsv")
anno_file <- file.path(
  nm2_root,
  "04_results", "figure5c_ncbi_clean", "figure5c_final", "figure5c_fresh_annotation_merged.csv"
)

out_dir <- file.path(nm2_root, "04_results", "figure5d_curated_examples_TABLE_FINAL_REVISED")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inst_file <- file.path(fig5b_dir, "instability_only_genes.tsv")
pool_file <- file.path(fig5b_dir, "pooled_only_genes.tsv")

for (fp in c(inst_file, pool_file, uv_file, anno_file)) {
  if (!file.exists(fp)) stop("Missing file: ", fp)
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clean_txt <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x[x %in% c("-", "/", "NA")] <- ""
  trimws(x)
}

shorten_annot <- function(x) {
  x <- clean_txt(x)
  x <- sub(";.*$", "", x)
  x <- sub(",.*$", "", x)
  x <- trimws(x)

  # manual cleanup for readability
  x[x == "Peptidase S1"] <- "Serine protease (trypsin family)"
  x[x == "P-loop containing nucleoside triphosphate hydrolase"] <- "Small GTPase (ARF/SAR type)"
  x
}

safe_col <- function(dt, candidates) {
  hit <- intersect(candidates, names(dt))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
inst <- fread(inst_file)
pool <- fread(pool_file)
uv   <- fread(uv_file)
anno <- fread(anno_file)

# ------------------------------------------------------------
# Standardize UV table
# ------------------------------------------------------------
if ("feature" %in% names(uv) && !("gene" %in% names(uv))) {
  setnames(uv, "feature", "gene")
}
if (!("gene" %in% names(uv))) stop("UV table has no gene column")
uv <- uv[!is.na(gene) & gene != ""]

# ------------------------------------------------------------
# Standardize gene-list tables
# ------------------------------------------------------------
for (obj_name in c("inst", "pool")) {
  obj <- get(obj_name)

  if (!("gene" %in% names(obj))) {
    cand <- safe_col(obj, c("feature", "Gene", "GENE"))
    if (is.na(cand)) stop(obj_name, " has no gene-like column")
    setnames(obj, cand, "gene")
  }

  obj <- obj[!is.na(gene) & gene != ""]
  assign(obj_name, obj)
}

# ------------------------------------------------------------
# Merge Figure 5B sets with UV metrics
# ------------------------------------------------------------
inst <- merge(inst, uv, by = "gene", all.x = TRUE)
pool <- merge(pool, uv, by = "gene", all.x = TRUE)

inst <- inst[!is.na(instability_score)]
pool <- pool[!is.na(delta_pool)]

# ------------------------------------------------------------
# Standardize annotation table
# ------------------------------------------------------------
if (!("gene" %in% names(anno))) {
  cand <- safe_col(anno, c("feature", "Gene", "GENE"))
  if (is.na(cand)) stop("Annotation table has no gene-like column")
  setnames(anno, cand, "gene")
}
anno <- anno[!is.na(gene) & gene != ""]

for (cc in intersect(c("ipr_desc", "ipr_go", "ipr_path", "pfam_hits"), names(anno))) {
  anno[[cc]] <- clean_txt(anno[[cc]])
}

anno_cols <- intersect(c("gene", "protein_id", "ipr_desc", "ipr_go", "ipr_path", "pfam_hits"), names(anno))
anno2 <- unique(anno[, ..anno_cols])

# ------------------------------------------------------------
# Merge annotation
# ------------------------------------------------------------
inst <- merge(inst, anno2, by = "gene", all.x = TRUE)
pool <- merge(pool, anno2, by = "gene", all.x = TRUE)

inst[, short_annot := fifelse(
  clean_txt(ipr_desc) != "", ipr_desc,
  fifelse(clean_txt(pfam_hits) != "", pfam_hits, "Uncharacterized protein")
)]
pool[, short_annot := fifelse(
  clean_txt(ipr_desc) != "", ipr_desc,
  fifelse(clean_txt(pfam_hits) != "", pfam_hits, "Uncharacterized protein")
)]

inst[, short_annot := shorten_annot(short_annot)]
pool[, short_annot := shorten_annot(short_annot)]

# ------------------------------------------------------------
# Tier 1: interpretable candidates
# ------------------------------------------------------------
inst_clean1 <- inst[clean_txt(ipr_desc) != "" | clean_txt(pfam_hits) != ""]
pool_clean1 <- pool[clean_txt(ipr_desc) != "" | clean_txt(pfam_hits) != ""]

inst_view1 <- inst_clean1[order(-instability_score), .(
  gene, instability_score, hetero_gap, sign_flip, delta_pool,
  short_annot, ipr_go, ipr_path, pfam_hits
)]

pool_view1 <- pool_clean1[order(-abs(delta_pool)), .(
  gene, delta_pool, instability_score, hetero_gap, sign_flip,
  short_annot, ipr_go, ipr_path, pfam_hits
)]

# ------------------------------------------------------------
# Tier 2 fallback: strong signal even if uncharacterized
# ------------------------------------------------------------
inst_fallback <- inst[
  short_annot == "Uncharacterized protein"
][
  order(-instability_score)
][
  sign_flip == TRUE | hetero_gap > quantile(hetero_gap, 0.9, na.rm = TRUE),
  .(gene, instability_score, hetero_gap, sign_flip, delta_pool,
    short_annot, ipr_go, ipr_path, pfam_hits)
]

pool_fallback <- pool[
  short_annot == "Uncharacterized protein"
][
  order(-abs(delta_pool))
][
  instability_score < quantile(instability_score, 0.5, na.rm = TRUE),
  .(gene, delta_pool, instability_score, hetero_gap, sign_flip,
    short_annot, ipr_go, ipr_path, pfam_hits)
]

inst_view2 <- unique(rbind(inst_view1, inst_fallback, fill = TRUE), by = "gene")
pool_view2 <- unique(rbind(pool_view1, pool_fallback, fill = TRUE), by = "gene")

inst_view2 <- inst_view2[!is.na(gene) & gene != ""]
pool_view2 <- pool_view2[!is.na(gene) & gene != ""]

# audit outputs
fwrite(inst_view2, file.path(out_dir, "inst_view2_ranked.tsv"), sep = "\t")
fwrite(pool_view2, file.path(out_dir, "pool_view2_ranked.tsv"), sep = "\t")

# ------------------------------------------------------------
# LOCKED final gene set
# ------------------------------------------------------------
inst_keep <- c("RvY_02420-1", "RvY_01282-1", "RvY_07877-1", "RvY_03909")
pool_keep <- c("RvY_03228-1", "RvY_06255-3")

inst_final <- inst_view2[gene %in% inst_keep]
pool_final <- pool_view2[gene %in% pool_keep]

if (nrow(inst_final) != length(inst_keep)) warning("Expected 4 instability genes, found ", nrow(inst_final))
if (nrow(pool_final) != length(pool_keep)) warning("Expected 2 pooled genes, found ", nrow(pool_final))

inst_final[, method := "Instability-only"]
pool_final[, method := "Pooled-only"]

# ------------------------------------------------------------
# Method-grounded behavior labels
# ------------------------------------------------------------
inst_final[, behavior := fifelse(
  sign_flip == TRUE,
  "Opposing regulation (sign flip)",
  "Magnitude divergence"
)]
pool_final[, behavior := "Consistent pooled response"]

# ------------------------------------------------------------
# Tight interpretations
# ------------------------------------------------------------
inst_final[, interpretation := fifelse(
  gene == "RvY_02420-1",
  "Opposing subgroup responses masked in pooled analysis",
  fifelse(
    gene == "RvY_01282-1",
    "Opposing subgroup behavior masked by pooling",
    fifelse(
      gene == "RvY_07877-1",
      "High instability independent of functional annotation",
      fifelse(
        gene == "RvY_03909",
        "High instability consistent with heterogeneous subgroup regulation",
        ""
      )
    )
  )
)]

pool_final[, interpretation := fifelse(
  gene == "RvY_03228-1",
  "Stable average signal prioritized by pooled ranking",
  fifelse(
    gene == "RvY_06255-3",
    "Consistent pooled signal with low internal conflict",
    ""
  )
)]

final_5c <- rbindlist(list(inst_final, pool_final), fill = TRUE)

# Clean display annotations
final_5c[, display_annot := short_annot]
final_5c[gene == "RvY_02420-1", display_annot := "Calycin-like protein"]
final_5c[gene == "RvY_01282-1", display_annot := "Serine protease (trypsin family)"]
final_5c[gene == "RvY_03228-1", display_annot := "Small GTPase (ARF/SAR type)"]

# reviewer-safe instability column
final_5c[, instability := sprintf("%.2f", instability_score)]

# Lock row order: annotated first within each method
final_5c[, method_order := fifelse(method == "Instability-only", 1L, 2L)]
final_5c[, annot_order := fifelse(display_annot == "Uncharacterized protein", 2L, 1L)]
final_5c[, abs_delta_pool := abs(delta_pool)]
setorder(final_5c, method_order, annot_order, -instability_score, -abs_delta_pool, gene)

final_5c_out <- final_5c[, .(
  gene,
  method,
  display_annot,
  behavior,
  instability,
  interpretation,
  instability_score,
  hetero_gap,
  sign_flip,
  delta_pool
)]

# ------------------------------------------------------------
# Export final curated table
# ------------------------------------------------------------
fwrite(final_5c_out, file.path(out_dir, "figure5d_final_curated_table.tsv"), sep = "\t")
write.csv(final_5c_out, file.path(out_dir, "figure5d_final_curated_table.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Build paper-style TABLE figure
# ------------------------------------------------------------
paper_tab <- copy(final_5c_out)[, .(
  Gene = gene,
  Method = method,
  Annotation = display_annot,
  Behavior = behavior,
  Instability = instability,
  Interpretation = interpretation
)]

# grouped fills
row_fills <- c("#F7F7F7", "#FFFFFF", "#F7F7F7", "#FFFFFF", "#EEF4FB", "#FFFFFF")

table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = 10, col = "black", hjust = 0, x = 0.02),
    bg_params = list(fill = row_fills[seq_len(nrow(paper_tab))], col = "grey80")
  ),
  colhead = list(
    fg_params = list(fontsize = 11, fontface = "bold", col = "black"),
    bg_params = list(fill = "#D9D9D9", col = "grey60")
  )
)

table_grob <- tableGrob(paper_tab, rows = NULL, theme = table_theme)

title_grob <- textGrob(
  "Figure 5D. Representative genes from instability-only and pooled-only rankings",
  gp = gpar(fontsize = 14, fontface = "bold")
)

subtitle_grob <- textGrob(
  "Instability-only genes capture subgroup-dependent or conflicting behavior, whereas pooled-only genes reflect stable average-prioritized signals.",
  gp = gpar(fontsize = 10)
)

foot_grob <- textGrob(
  "Instability-only examples include annotated candidates and highly unstable uncharacterized proteins, indicating masked structured biology beyond stable pooled responses.",
  gp = gpar(fontsize = 9)
)

# PNG
png(file.path(out_dir, "Figure5d_curated_examples_TABLE_FINAL_REVISED.png"),
    width = 3000, height = 800, res = 260)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  4, 1,
  heights = unit(c(0.09, 0.07, 0.76, 0.08), "npc")
)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(title_grob)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(subtitle_grob)
upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
grid.draw(table_grob)
upViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
grid.draw(foot_grob)
upViewport(2)
dev.off()

# PDF
pdf(file.path(out_dir, "Figure5d_curated_examples_TABLE_FINAL_REVISED.pdf"),
    width = 14.2, height = 3)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  4, 1,
  heights = unit(c(0.09, 0.07, 0.76, 0.08), "npc")
)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(title_grob)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(subtitle_grob)
upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
grid.draw(table_grob)
upViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
grid.draw(foot_grob)
upViewport(2)
dev.off()

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------
cat("Figure 5d TABLE final revised panel generated.\n\n")
cat("Output directory:\n", out_dir, "\n\n", sep = "")
cat("Final selected genes:\n")
print(final_5c_out)
