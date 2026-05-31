# Repository Map

This document provides a direct mapping between manuscript figures, analysis scripts, generated outputs, and supplementary materials.

---

# Main Figures

## Figure 1

**Title:** Pooling masks structured regulatory signals, whereas instability reveals them.

**Output**

```text
03_figures/main/Fig1.tif
```

**Description**
Conceptual illustration of pooling-induced signal attenuation and instability-based recovery.

---

## Figure 2

### Figure 2a

**Synthetic expression regimes**

Script:

```text
02_methods/R/core_pipeline/figure2a_synthetic_expression_regimes.R
```

### Figure 2b

**Pooled failure under subgroup conflict**

Script:

```text
02_methods/R/figure_generation/figure2_pooled_failure.R
```

### Figure 2c

**Instability enrichment**

Script:

```text
02_methods/R/figure_generation/figure2_instability_enrichment.R
```

Output:

```text
03_figures/main/Fig2.tiff
```

---

## Figure 3

### Figure 3a

Script:

```text
02_methods/R/figure_generation/figure3_uv_quadrant.R
```

### Figure 3b

Script:

```text
02_methods/R/figure_generation/figure3_instability_vs_variance.R
```

### Figure 3c

Script:

```text
02_methods/R/figure_generation/figure3_instability_vs_disagreement.R
```

### Figure 3d

Script:

```text
02_methods/R/figure_generation/figure3_gene_examples.R
```

Output:

```text
03_figures/main/Fig3.tiff
```

---

## Figure 4

### Figure 4a

Script:

```text
02_methods/R/figure_generation/figure4_des_support.R
```

### Figure 4b

Script:

```text
02_methods/R/figure_generation/figure4_gam_temporal_extension.R
```

### Figure 4c-d

Script:

```text
02_methods/R/figure_generation/figure4_cross_condition_panel.R
```

Output:

```text
03_figures/main/Fig4.png
```

---

## Figure 5

### Figure 5a

Script:

```text
02_methods/R/figure_generation/figure5_rank_displacement.R
```

### Figure 5b

Script:

```text
02_methods/R/figure_generation/figure5_topN_overlap.R
```

### Figure 5c-e

Generated through annotation workflow:

```text
02_methods/R/annotation_pipeline/
```

Scripts:

```text
00_validate_mapping_integrity.R
01_build_ncbi_reference_map_for_figure5c.R
02_match_figure5c_genes_to_ncbi_reference.R
03_extract_matched_ncbi_proteins_for_figure5c.R
04_merge_annotation_and_make_figure5c.R
05_term_enrichment_for_figure5c.R
```

Output:

```text
03_figures/main/Fig5.tif
```

---

## Figure 6

Script:

```text
02_methods/R/figure_generation/figure6_variance_disentanglement.R
```

Output:

```text
03_figures/main/Fig6.tif
```

---

# Supplementary Figures

## Figure S1

Output:

```text
03_figures/supplementary/FigS1.png
```

Description:
DES magnitude-divergence support.

---

## Figure S2

Output:

```text
03_figures/supplementary/FigS2.png
```

Description:
GAM temporal-divergence support.

---

# Supplementary Tables

## Table S1

```text
05_supplementary/supplementary_tables/TableS1_metadata.csv
```

Sample metadata and experimental structure.

---

## Table S2

```text
05_supplementary/supplementary_tables/TableS2_representative_masked_genes.csv
```

Representative masked genes showing subgroup effect cancellation.

---

## Table S3

```text
05_supplementary/supplementary_tables/TableS3_functional_annotations.xlsx
```

Functional annotation of instability-recovered masked genes.

---

# Instability vs Variance Package

Location:

```text
05_supplementary/instability_vs_variance/
```

Contents:

* Correlation analyses
* Regression outputs
* Synthetic benchmarks
* Supporting raw tables

Purpose:

Provide supporting evidence that instability is not reducible to variance alone and address potential reviewer questions regarding independence from variance.

---

# Development History

Historical exploratory scripts are retained for transparency:

```text
archive/
```

These scripts were not required to reproduce the final manuscript figures and results.

