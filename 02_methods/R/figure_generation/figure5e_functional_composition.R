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

# ============================================================
# Paths
# ============================================================
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
pfam_file <- file.path(ann_dir, "figure5e_pfam.domtblout")
ipr_file  <- file.path(ann_dir, "figure5e_interproscan.tsv")

out_dir <- file.path(nm2_root, "04_results", "figure5e_FINAL")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# outputs
masked_out        <- file.path(out_dir, "figure5e_masked_genes.tsv")
mapping_out       <- file.path(out_dir, "figure5e_gene_protein_mapping.tsv")
annot_table_out   <- file.path(out_dir, "figure5e_raw_annotation_table.tsv")
counts_out        <- file.path(out_dir, "figure5e_category_counts.tsv")
coverage_out      <- file.path(out_dir, "figure5e_annotation_coverage.tsv")
expert_override_out <- file.path(out_dir, "figure5e_expert_manual_reclassification.tsv")
expert_reclassified_out <- file.path(out_dir, "figure5e_expert_reclassified_genes.tsv")
complex_remaining_out <- file.path(out_dir, "figure5e_remaining_complex_unresolved.tsv")
plot_png          <- file.path(out_dir, "figure5e_FINAL.png")
plot_pdf          <- file.path(out_dir, "figure5e_FINAL.pdf")

for (fp in c(uv_file, cds_fna, pfam_file, ipr_file)) {
  if (!file.exists(fp)) stop("Missing file: ", fp)
}

# ============================================================
# Helpers
# ============================================================
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

assign_primary_category <- function(txt) {
  txt <- tolower(clean_txt(txt))
  txt <- gsub("\\s+", " ", txt)

  if (txt == "") return("Uncharacterized")

  if (grepl("receptor|gpcr|kinase|phosphatase|g-protein|signal|pdz|lrr|leucine rich repeat|hormone|neuropeptide|peptide", txt)) {
    return("Signaling/Receptor")
  }
  if (grepl("protease|peptidase|trypsin|serine protease|cathepsin|metalloprotease", txt)) {
    return("Proteolysis")
  }
  if (grepl("transporter|channel|carrier|membrane|slc|sulp|stas|mucin", txt)) {
    return("Transport/Membrane")
  }
  if (grepl("transcription|dna-binding|rna-binding|regulator|chromatin|zinc finger|polymerase|cdc27|cdc45|endonuclease", txt)) {
    return("Regulation/Nucleic-acid associated")
  }
  if (grepl("dehydrogenase|synthetase|hydrolase|transferase|oxidoreductase|enzyme|metabol|esterase|epimerase|phosphatase", txt)) {
    return("Metabolism/Enzyme")
  }
  if (grepl("actin|tubulin|cytoskeleton|structural|motor protein|coiled-coil|capsid|nucleocapsid|phage|viral", txt)) {
    return("Structural/Other")
  }
  if (grepl("ribosome|ribosomal|translation|splice|elongation|rna processing", txt)) {
    return("RNA / translation")
  }
  if (grepl("duf|c2orf69", txt)) {
    return("DUF / unknown conserved")
  }
  if (grepl("binding|conotoxin|self-incompatibility", txt)) {
    return("Binding proteins")
  }
  if (grepl("repeat|scaffold|ankyrin|wd40|armadillo|tpr", txt)) {
    return("Repeat / scaffold proteins")
  }
  if (grepl("domain|superfamily|family|motif|fold", txt)) {
    return("Generic domain proteins")
  }

  return("Complex proteins, unresolved class")
}

collapse_for_figure <- function(cat) {
  case_when(
    cat %in% c(
      "Complex proteins, unresolved class",
      "Complex domain proteins",
      "Conserved complex proteins",
      "Unresolved: complex unknown",
      "Unresolved: structural/disordered",
      "Unresolved: conserved generic domain",
      "Unresolved: catalytic-like",
      "Unresolved: binding/regulatory-like",
      "Unresolved: membrane-associated",
      "Unresolved: repeat/scaffold-rich",
      "Unresolved: no detectable feature"
    ) ~ "Complex proteins, unresolved class",
    TRUE ~ cat
  )
}

# ============================================================
# Expert manual override table
# ============================================================
expert_override <- tibble::tribble(
  ~pattern, ~expert_group, ~expert_note,
  "pro-corazonin", "Signaling/Receptor", "Neuropeptide/hormone precursor; stress-physiology signaling component",
  "myotubularin", "Signaling/Receptor", "Signaling-associated phosphatase; enzyme-related but retained under signaling",
  "alpha conotoxin", "Signaling/Receptor", "Bioactive binding peptide/toxin-like precursor",
  "ribosomal protein s2", "RNA / translation", "Core ribosomal translation component",
  "dna/rna non-specific endonuclease", "Regulation/Nucleic-acid associated", "DNA/RNA nuclease-associated feature",
  "endonuclease", "Regulation/Nucleic-acid associated", "DNA/RNA nuclease-associated feature",
  "arthropod cardioacceleratory peptide", "Signaling/Receptor", "Neuropeptide/hormonal signaling peptide",
  "dna polymerase subunit cdc27", "Regulation/Nucleic-acid associated", "DNA polymerase-associated regulatory feature",
  "cdc27", "Regulation/Nucleic-acid associated", "Cell-cycle/DNA-associated regulatory component",
  "putative sensor; cdc45", "Regulation/Nucleic-acid associated", "DNA replication-associated sensor",
  "cdc45", "Regulation/Nucleic-acid associated", "DNA replication-associated sensor",
  "mucin-like glycoprotein", "Structural/Other", "Extracellular/surface protective protein; possible membrane/surface role",
  "phage pam3 gp32-like protein", "Structural/Other", "Viral-like structural annotation",
  "pam3 gp32", "Structural/Other", "Viral-like structural annotation",
  "coronavirus nucleocapsid", "Structural/Other", "Viral-like structural annotation",
  "plant self-incompatibility protein s1", "Binding proteins", "Plant-like binding annotation",
  "pigment-dispersing hormone", "Signaling/Receptor", "Neuropeptide hormone; circadian/behavioral signaling component",
  "pdh", "Signaling/Receptor", "Pigment-dispersing hormone-like neuropeptide annotation",
  "pectinacetylesterase", "Metabolism/Enzyme", "Possible microbiome/contamination-associated enzyme",
  "c2orf69", "DUF / unknown conserved", "Conserved unknown protein",
  "phi6 bacteriophage p8 capsid protein", "Structural/Other", "Viral-like structural annotation",
  "phi6", "Structural/Other", "Viral-like structural annotation",
  "diaminopimelate epimerase", "Metabolism/Enzyme", "Possible bacterial contamination-associated enzyme"
)

write_tsv(expert_override, expert_override_out)

apply_expert_override <- function(raw_annotation, current_group) {
  txt <- tolower(clean_txt(raw_annotation))

  hit <- expert_override %>%
    filter(str_detect(txt, fixed(tolower(pattern))))

  if (nrow(hit) == 0) return(current_group)
  hit$expert_group[1]
}

apply_expert_note <- function(raw_annotation) {
  txt <- tolower(clean_txt(raw_annotation))

  hit <- expert_override %>%
    filter(str_detect(txt, fixed(tolower(pattern))))

  if (nrow(hit) == 0) return("")
  hit$expert_note[1]
}

# ============================================================
# 1) Define masked genes from UV table
# ============================================================
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
  filter(
    abs(delta_pool) <= q25_pool,
    instability_score >= q90_inst
  ) %>%
  distinct(gene, .keep_all = TRUE)

write_tsv(masked, masked_out)
cat("Masked genes:", nrow(masked), "\n")

# ============================================================
# 2) Parse CDS headers for gene -> protein mapping
# ============================================================
cat("Parsing cds_from_genomic.fna headers...\n")
cds_lines <- readLines(cds_fna, warn = FALSE)
cds_headers <- cds_lines[startsWith(cds_lines, ">")]

cds_map <- tibble(
  raw_header = cds_headers,
  gene = extract_tag(cds_headers, "gene"),
  locus_tag = extract_tag(cds_headers, "locus_tag"),
  ncbi_protein_desc = extract_tag(cds_headers, "protein"),
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

cat("CDS entries with protein_id:", nrow(cds_map), "\n")

# ============================================================
# 3) Map masked genes to proteins
# ============================================================
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
  select(
    gene,
    gene_base,
    protein_id,
    ncbi_protein_desc,
    match_type,
    delta_pool,
    instability_score
  )

write_tsv(mapped, mapping_out)
cat("Mapped masked genes:", nrow(mapped), "\n")

# ============================================================
# 4) Parse Pfam
# ============================================================
cat("Parsing Pfam output...\n")
pfam_lines <- readLines(pfam_file, warn = FALSE)
pfam_lines <- pfam_lines[!startsWith(pfam_lines, "#")]

pfam_df <- lapply(pfam_lines, function(ln) {
  z <- strsplit(trimws(ln), "\\s+")[[1]]
  if (length(z) < 23) return(NULL)

  data.frame(
    protein_id = z[4],
    pfam_name  = z[1],
    pfam_acc   = z[2],
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
    pfam_hit_count = n(),
    .groups = "drop"
  )

cat("Pfam-annotated proteins:", nrow(pfam_sum), "\n")

# ============================================================
# 5) Parse InterProScan
# ============================================================
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
    interpro_analysis = collapse_unique(analysis),
    interpro_sig_acc  = collapse_unique(sig_acc),
    interpro_sig_desc = collapse_unique(sig_desc),
    interpro_acc      = collapse_unique(ipr_acc),
    interpro_desc     = collapse_unique(ipr_desc),
    interpro_hit_count = n(),
    .groups = "drop"
  )

cat("InterPro-annotated proteins:", nrow(ipr_sum), "\n")

# ============================================================
# 6) Merge, automatically classify, then expert-reclassify
# ============================================================
cat("Merging annotation evidence and applying expert curation...\n")

annot_tbl <- mapped %>%
  left_join(pfam_sum, by = "protein_id") %>%
  left_join(ipr_sum, by = "protein_id") %>%
  rowwise() %>%
  mutate(
    raw_annotation = collapse_unique(c(
      ncbi_protein_desc,
      pfam_desc,
      interpro_sig_desc,
      interpro_desc,
      pfam_name,
      interpro_analysis
    )),
    has_annotation = raw_annotation != "",

    broad_group_auto = ifelse(
      has_annotation,
      assign_primary_category(raw_annotation),
      "Uncharacterized"
    ),

    expert_reclassified_group = apply_expert_override(
      raw_annotation,
      broad_group_auto
    ),

    expert_reclassification_note = apply_expert_note(raw_annotation),

    manual_reclassified = expert_reclassification_note != "",

    broad_group_precollapse = expert_reclassified_group,

    broad_group_figure = collapse_for_figure(broad_group_precollapse)
  ) %>%
  ungroup() %>%
  arrange(desc(instability_score), gene)

write_tsv(annot_tbl, annot_table_out)

expert_reclassified_tbl <- annot_tbl %>%
  filter(manual_reclassified) %>%
  select(
    gene,
    protein_id,
    raw_annotation,
    broad_group_auto,
    broad_group_figure,
    expert_reclassification_note,
    delta_pool,
    instability_score
  )

write_tsv(expert_reclassified_tbl, expert_reclassified_out)

remaining_complex_tbl <- annot_tbl %>%
  filter(broad_group_figure == "Complex proteins, unresolved class") %>%
  select(
    gene,
    protein_id,
    raw_annotation,
    broad_group_auto,
    broad_group_figure,
    delta_pool,
    instability_score
  )

write_tsv(remaining_complex_tbl, complex_remaining_out)

cat("Expert-reclassified entries:", nrow(expert_reclassified_tbl), "\n")
cat("Remaining complex unresolved entries:", nrow(remaining_complex_tbl), "\n")

# ============================================================
# 7) Coverage
# ============================================================
coverage <- annot_tbl %>%
  summarise(
    total = n(),
    annotated = sum(has_annotation),
    unannotated = sum(!has_annotation),
    annotated_pct = round(100 * annotated / total, 1),
    expert_reclassified = sum(manual_reclassified),
    remaining_complex_unresolved = sum(broad_group_figure == "Complex proteins, unresolved class")
  )

write_tsv(coverage, coverage_out)

cat("\nAnnotation coverage:\n")
print(coverage)

# ============================================================
# 8) Final counts for figure
# ============================================================
counts <- annot_tbl %>%
  count(broad_group_figure, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1))

write_tsv(counts, counts_out)

cat("\nFinal figure category counts:\n")
print(counts)

# ============================================================
# 9) Final clean plot
# ============================================================
category_colors <- c(
  "Signaling/Receptor" = "#D55E00",
  "Metabolism/Enzyme" = "#E69F00",
  "Transport/Membrane" = "#56B4E9",
  "Regulation/Nucleic-acid associated" = "#009E73",
  "Proteolysis" = "#CC79A7",
  "Binding proteins" = "#0072B2",
  "Structural/Other" = "#999999",
  "RNA / translation" = "#F0E442",
  "DUF / unknown conserved" = "#BDBDBD",
  "Generic domain proteins" = "#8DD3C7",
  "Repeat / scaffold proteins" = "#80B1D3",
  "Complex proteins, unresolved class" = "#636363",
  "Uncharacterized" = "#CCCCCC"
)

counts$color_group <- ifelse(
  counts$broad_group_figure %in% names(category_colors),
  counts$broad_group_figure,
  "Uncharacterized"
)

p <- ggplot(
  counts,
  aes(
    x = reorder(broad_group_figure, n),
    y = n,
    fill = color_group
  )
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = category_colors, drop = FALSE) +
  coord_flip(clip = "off") +
  labs(
    title = "Functional composition of instability-recovered masked genes",
    x = "Broad functional class",
    y = "Number of genes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13,
      margin = margin(b = 12)
    ),
    plot.title.position = "plot",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 16, r = 12, b = 10, l = 12)
  )

try(
  ggsave(plot_png, p, width = 7.2, height = 5.4, dpi = 600, bg = "white", device = "png"),
  silent = TRUE
)

ggsave(plot_pdf, p, width = 7.2, height = 5.4, bg = "white", device = "pdf")

if (!file.exists(plot_pdf)) {
  stop("PDF export failed for Figure 5e: ", plot_pdf)
}

cat("\nFinal Figure 5E outputs written to:\n")
cat(out_dir, "\n")
cat("Main figure PNG:", plot_png, "\n")
cat("Main figure PDF:", plot_pdf, "\n")
cat("Raw annotation table:", annot_table_out, "\n")
cat("Expert override table:", expert_override_out, "\n")
cat("Expert-reclassified genes:", expert_reclassified_out, "\n")
cat("Remaining unresolved complex genes:", complex_remaining_out, "\n")
cat("Category counts:", counts_out, "\n")
cat("Coverage summary:", coverage_out, "\n")
