# Instability vs Variance

## Purpose

A central claim of the NM2 framework is that instability captures structured disagreement among partitions and is not reducible to variance alone.

Because variance is commonly used to quantify dispersion in transcriptomic datasets, an important methodological question is whether instability simply represents a transformed variance measure. The analyses in this directory were assembled to evaluate that possibility directly.

This package provides supporting evidence that instability and variance capture related but distinct properties of gene expression behavior.

---

## Question Addressed

> Is instability simply another measure of variance?

The analyses in this directory evaluate this question using both real and synthetic datasets.

---

## Directory Contents

```text
instability_vs_variance/

README.md

correlation_analysis/
regression_outputs/
synthetic_benchmarks/
raw_tables/
```

---

## correlation_analysis

Contains exploratory analyses evaluating the relationship between instability and variance.

Typical contents:

```text
correlation_analysis/

instability_vs_variance_scatter.png
correlation_summary.csv
```

These analyses assess whether instability exhibits a simple linear or monotonic relationship with variance.

---

## regression_outputs

Contains statistical models evaluating the contribution of instability and variance to masking behavior.

Typical contents:

```text
regression_outputs/

logistic_regression_results.csv
coefficient_table.csv
model_summary.txt
```

These analyses correspond to the regression framework described in the manuscript and evaluate whether instability contributes information beyond variance alone.

---

## synthetic_benchmarks

Contains synthetic validation datasets and summary outputs.

Typical contents:

```text
synthetic_benchmarks/

coherent/
phase/
opposing/

synthetic_summary.csv
```

These benchmarks evaluate how instability and variance behave under controlled response regimes.

The goal is to determine whether the two quantities scale proportionally across scenarios or capture different structural properties.

---

## raw_tables

Contains source tables used for instability-versus-variance analyses.

Typical contents:

```text
raw_tables/

instability_metrics.tsv
variance_metrics.tsv
masked_gene_table.tsv
```

These files provide direct access to the underlying values used in the analyses.

---

## Interpretation

Variance quantifies dispersion.

Instability quantifies disagreement among structured partitions.

Although both quantities may increase under some forms of heterogeneous regulation, they are not expected to behave identically across all scenarios.

The analyses contained in this directory were generated to support transparent evaluation of the relationship between these measures and to facilitate independent verification of the results reported in the manuscript.

---

## Related Manuscript Sections

Main manuscript:

* Figure 3
* Figure 6
* Discussion

Supplementary analyses contained here provide additional supporting evidence beyond the figures presented in the manuscript.

