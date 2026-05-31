# NM2: Instability-Based Detection of Structured Heterogeneity in Transcriptomic Data

## Overview

NM2 (Instability Signal) is a computational framework designed to identify structured heterogeneity that can be attenuated or underrepresented by pooled transcriptomic analyses. Rather than focusing solely on average effect size, NM2 quantifies disagreement among biologically meaningful partitions such as subgroups, exposure regimes, or timepoints.

The framework was developed and evaluated using RNA-seq datasets from *Ramazzottius varieornatus* under multiple stress conditions, including ultraviolet radiation (UV), desiccation (DES), gamma radiation (GAM), and low temperature (LT).

The accompanying manuscript demonstrates that instability identifies genes exhibiting structured subgroup-dependent or temporally heterogeneous responses that may be attenuated under pooled analysis.

---

## Repository Structure

```text
NM2_instability_signal/
├── 01_data/
├── 02_methods/
├── 03_figures/
├── 04_results/
├── 05_supplementary/
├── manuscript/
└── docs/
```

---

## Core Components

### `02_methods/R/core_pipeline/`

Core analytical workflow for pooled and structured analyses, instability calculation, heterogeneity assessment, and summary table generation.

### `02_methods/R/figure_generation/`

Scripts used to generate manuscript figures.

### `02_methods/R/figure5c_annotation_workflow/`

Supporting workflow used for biological annotation and functional characterization of instability-prioritized genes.

### `02_methods/R/archive/`

Development and exploratory scripts retained for transparency but not required for manuscript reproduction.

---

## Main Figures

| Figure | Description                                             |
| ------ | ------------------------------------------------------- |
| Fig. 1 | Conceptual illustration of pooling-induced masking      |
| Fig. 2 | Synthetic validation of instability                     |
| Fig. 3 | UV instability and masked regulatory signals            |
| Fig. 4 | Cross-condition heterogeneity patterns                  |
| Fig. 5 | Prioritization divergence and biological interpretation |
| Fig. 6 | Instability versus variance disentanglement             |

---

## Reproducing the Manuscript

The complete manuscript workflow can be reproduced from a clean repository checkout using:

```bash
Rscript 02_methods/R/run_nm2_pipeline.R
```

All manuscript figures, tables, and intermediate outputs are written to:

```text
04_results/
```

The master runner executes analytical scripts as isolated `Rscript --vanilla` processes.

Resources requiring external software are supplied as bundled precomputed files to ensure reproducibility.

### Figure 5e Annotation Resources

Figure 5e uses supplied InterProScan and Pfam/HMMER outputs included in:

```text
01_data/precomputed/figure5e_annotation/
```

The manuscript pipeline does not rerun external annotation software.

### Figure 6 Variance Package

Figure 6 source files are bundled as precomputed resources in:

```text
01_data/precomputed/figure6/
```

---

## Supplementary Materials

Supplementary figures, tables, and supporting analyses are available in:

```text
05_supplementary/
```

### Instability versus Variance Package

```text
05_supplementary/instability_vs_variance/
```

This directory contains:

* Correlation analyses
* Regression models
* Synthetic benchmarks
* Supporting raw tables

---

## Manuscript

The manuscript associated with this repository is available in:

```text
manuscript/
```

---

## Citation

If you use this repository, please cite the associated publication.

Citation information will be updated upon acceptance or publication of the associated manuscript.

