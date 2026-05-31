# Repository Structure

This document describes the organization of the NM2 Instability Signal repository and the purpose of each directory.

---

# Top-Level Layout

```text
NM2_instability_signal/

README.md
REPOSITORY_MAP.md

01_data/
02_methods/
03_results/
04_figures/
05_supplementary/

manuscript/
docs/
tests/
```

---

# 01_data

Contains processed input datasets and metadata used for analysis.

```text
01_data/
├── metadata/
├── processed/
├── synthetic/
└── examples/
```

### metadata

Sample metadata tables and condition definitions.

### processed

Processed expression matrices and derived datasets used by the analysis workflow.

### synthetic

Synthetic benchmark datasets used to evaluate instability behavior under controlled scenarios.

### examples

Small example datasets for workflow demonstration and testing.

---

# 02_methods

Contains all analysis scripts and computational workflows.

```text
02_methods/
└── R/
    ├── core_pipeline/
    ├── figure_generation/
    ├── annotation_pipeline/
    └── run_nm2_pipeline.R
```

### core_pipeline

Core analytical workflow including pooled analysis, structured analysis, instability calculation, masking analysis, and summary statistics.

### figure_generation

Scripts used to generate manuscript figures.

### annotation_pipeline

Functional annotation workflows used for Figure 5 and Supplementary Table S3.

### run_nm2_pipeline.R

Master script for reproducing major analyses and figures.

---

# 03_results

Contains generated outputs produced by analysis scripts.

```text
03_results/
├── figures/
├── statistics/
└── tables/
```

### figures

Intermediate figure outputs.

### statistics

Statistical summaries and model outputs.

### tables

Generated result tables used in the manuscript.

---

# 04_figures

Contains publication-ready figures.

```text
04_figures/
├── main/
└── supplementary/
```

### main

Figures included in the main manuscript.

### supplementary

Figures included in supplementary materials.

---

# 05_supplementary

Contains supplementary datasets, supporting analyses, and reviewer-support materials.

```text
05_supplementary/
├── supplementary_figures/
├── supplementary_tables/
└── instability_vs_variance/
```

### supplementary_figures

Supplementary figures and captions.

### supplementary_tables

Supplementary tables and associated metadata.

### instability_vs_variance

Supporting analyses evaluating the relationship between instability and variance, including correlation analyses, regression models, synthetic benchmarks, and raw data tables.

This section is intended to support reproducibility and address methodological questions regarding the distinction between instability and variance.

---

# manuscript

Contains manuscript files and publication materials.

```text
manuscript/
├── main_text/
├── supplementary/
└── responses/
```

### main_text

Submitted manuscript versions.

### supplementary

Supplementary information accompanying the manuscript.

### responses

Reviewer responses and revision materials.

---

# docs

Repository documentation.

```text
docs/
├── repository_structure.md
└── ...
```

Includes repository organization, workflow descriptions, and reproducibility notes.

---

# tests

Contains validation and reproducibility checks.

```text
tests/
├── figure_reproducibility/
└── script_validation/
```

These tests are intended to verify that the analysis pipeline executes successfully and reproduces expected outputs.

---

# Archive Scripts

Exploratory and development scripts not required for manuscript reproduction are retained separately in archive directories for transparency.

These files are not part of the minimal reproducible workflow and are not required to regenerate manuscript results.

