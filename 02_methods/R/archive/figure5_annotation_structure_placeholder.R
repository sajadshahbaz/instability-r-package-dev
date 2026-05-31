#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ==============================
# INPUT FILES (EDIT THESE PATHS)
# ==============================

# Your main results table (must include: gene_id, delta_pool, instability)
INPUT_RESULTS <- "results_table.tsv"

# Your merged annotation table from previous pipeline
INPUT_ANNOT   <- "annotation_table.tsv"

# Output
OUT_FIG <- "Figure5E_functional_categories.png"
OUT_TAB <- "Figure5E_category_counts.tsv"

# ==============================
# LOAD DATA
# ==============================

cat("Loading data...\n")

df <- fread(INPUT_RESULTS)
annot <- fread(INPUT_ANNOT)

# sanity check
required_cols <- c("gene_id", "delta_pool", "instability")
if (!all(required_cols %in% colnames(df))) {
  stop("Missing required columns in results table")
}

# ==============================
# STEP 1: DEFINE MASKED GENES
# ==============================

cat("Selecting masked genes...\n")

Q25_pool <- quantile(abs(df$delta_pool), 0.25, na.rm = TRUE)
Q90_inst <- quantile(df$instability, 0.90, na.rm = TRUE)

masked <- df[
  abs(delta_pool) <= Q25_pool &
  instability >= Q90_inst
]

cat("Masked genes:", nrow(masked), "\n")

# ==============================
# STEP 2: MERGE ANNOTATION
# ==============================

cat("Merging annotation...\n")

masked_annot <- merge(
  masked,
  annot,
  by = "gene_id",
  all.x = TRUE
)

# ==============================
# STEP 3: CATEGORY ASSIGNMENT
# ==============================

assign_category <- function(desc, pfam, interpro, eggnog) {

  text <- tolower(paste(desc, pfam, interpro, eggnog, collapse = " "))

  if (grepl("receptor|gpcr|membrane", text)) return("Receptor/Membrane")
  if (grepl("signal|kinase|phosphat", text)) return("Signaling")
  if (grepl("stress|heat|shock|oxidative", text)) return("Stress response")
  if (grepl("transport|channel|carrier", text)) return("Transport")
  if (grepl("protease|ubiquitin|degrad", text)) return("Proteolysis")
  if (grepl("transcription|rna|dna binding", text)) return("Regulation")
  if (grepl("metabolic|enzyme", text)) return("Metabolism")
  if (grepl("repair|dna repair", text)) return("DNA repair")
  if (grepl("cytoskeleton|structur", text)) return("Structural")

  return("Uncharacterized")
}

cat("Assigning functional categories...\n")

masked_annot[, category := mapply(
  assign_category,
  description,
  pfam,
  interpro,
  eggnog
)]

# ==============================
# STEP 4: COUNT CATEGORIES
# ==============================

cat("Counting categories...\n")

cat_counts <- masked_annot[
  , .N, by = category
][order(-N)]

cat_counts[, percent := round(100 * N / sum(N), 2)]

# save table
fwrite(cat_counts, OUT_TAB, sep = "\t")

# ==============================
# STEP 5: PLOT FIGURE 5E
# ==============================

cat("Plotting figure...\n")

p <- ggplot(cat_counts, aes(x = reorder(category, N), y = N)) +
  geom_bar(stat = "identity", fill = "#2C7BB6") +
  coord_flip() +
  labs(
    x = "Functional category",
    y = "Number of genes",
    title = "Functional composition of instability-recovered masked genes"
  ) +
  theme_minimal(base_size = 12)

ggsave(OUT_FIG, p, width = 6, height = 4, dpi = 300)

cat("Done.\n")
