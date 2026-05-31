#!/usr/bin/env bash
set -euo pipefail

test -f 03_figures/main/Fig1_conceptual_masking.tiff
test -f 03_figures/main/Fig2_synthetic_validation.tiff
test -f 03_figures/main/Fig3_uv_instability.tiff
test -f 03_figures/main/Fig4_cross_condition_specificity.tiff
test -f 03_figures/main/Fig5_prioritization_divergence.tiff
test -f 03_figures/main/Fig6_variance_disentanglement.tiff

test -f 04_results/figure5a_rank_displacement/figure5a_rank_displacement_locked.pdf
test -f 04_results/figure5b_topN_overlap/figure5b_topN_overlap_locked.pdf
test -f 04_results/figure5c_selection_stability_FINAL/Figure5C_selection_stability_FINAL.pdf
test -f 04_results/figure5d_curated_examples_TABLE_FINAL_REVISED/Figure5d_curated_examples_TABLE_FINAL_REVISED.pdf
test -f 04_results/figure5e_FINAL/figure5e_FINAL.pdf

test -f 04_results/figure5e_global_functional_structure/annotation_from_scratch/figure5e_pfam.domtblout
test -f 04_results/figure5e_global_functional_structure/annotation_from_scratch/figure5e_interproscan.tsv

echo "Manuscript output audit passed."
