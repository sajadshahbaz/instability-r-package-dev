# NM2: Instability-Based Detection of Structured Heterogeneity in Transcriptomic Data

## Overview

NM2 (Instability Signal) is a computational framework designed to identify structured heterogeneity that can be attenuated or underrepresented by pooled transcriptomic analyses. Rather than focusing solely on average effect size, NM2 quantifies disagreement among biologically meaningful partitions such as subgroups, exposure regimes, or timepoints.

The framework was developed and evaluated using RNA-seq datasets from the tardigrade *Ramazzottius varieornatus* under multiple stress conditions, including ultraviolet radiation (UV), desiccation (DES), gamma radiation (GAM), and low temperature (LT).

The accompanying manuscript demonstrates that instability identifies genes exhibiting structured subgroup-dependent or temporally heterogeneous responses that may be attenuated under pooled analysis.

---

## Repository Structure

```text
NM2_instability_signal/

01_data/
02_methods/
03_figures/
04_results/
05_supplementary/
manuscript/
docs/
```

### Core Components

**02_methods/R/core_pipeline/**
Core analytical workflow for pooled and structured analyses, instability calculation, heterogeneity assessment, and summary table generation.

**02_methods/R/figure_generation/**
Scripts used to generate manuscript figures.

**02_methods/R/annotation_pipeline/**
Workflow used for annotation and functional characterization of instability-prioritized genes.

**archive/**
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

### Step 1

Run the core workflow:

```r
source("02_methods/R/core_pipeline/...")
```

### Step 2

Generate manuscript figures:

```r
source("02_methods/R/figure_generation/...")
```

### Step 3

Generate annotation outputs:

```r
source("02_methods/R/annotation_pipeline/...")
```

Outputs will be written to:

```text
03_figures/
04_results/
```

---

## Supplementary Materials

Supplementary figures, tables, and supporting analyses are available in:

```text
05_supplementary/
```

Of particular interest:

```text
05_supplementary/instability_vs_variance/
```

which contains analyses evaluating the relationship between instability and variance, including:

* correlation analyses
* regression models
* synthetic benchmarks
* supporting raw tables

---

## Manuscript

The manuscript associated with this repository is located in:

```text
manuscript/
```

---

## Citation

If you use this repository, please cite the associated publication.

Citation information will be updated following publication.

