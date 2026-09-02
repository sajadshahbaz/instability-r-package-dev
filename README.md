# instability

`instability` provides the validated ultraviolet (UV) and desiccation (DES)
Instability computation for already-preprocessed variance-stabilized expression
data and associated metadata.

For each condition, the package compares treatment-minus-control effects
calculated after pooling samples with effects calculated separately within two
subgroups. It reports feature-level heterogeneity, disagreement between pooled
and stratified effects, ranking changes, and condition-level summaries.

The implementation intentionally preserves the scientific behavior of the
authoritative UV/DES implementation. It is not a general framework for
arbitrary condition names or subgroup designs.

## Installation

Version 0.1.1 requires R 4.5.1 or later and imports `data.table` and `stats`.

Install a downloaded v0.1.1 source-package tarball:

```r
install.packages(
  "instability_0.1.1.tar.gz",
  repos = NULL,
  type = "source"
)
```

Or install from a local source checkout:

```text
R CMD INSTALL /path/to/instability-r-package-dev
```

The release tarball must be obtained from the package repository's v0.1.1
release. The package is not claimed to be available from CRAN.

Problems may be reported at:

https://github.com/sajadshahbaz/instability-r-package-dev/issues

## Public API

The package exports exactly one function:

```r
compute_uv_des_instability(vst, metadata)
```

There are no public scientific tuning arguments in version 0.1.1.

## VST input

`vst` must be a `data.frame` or `data.table`.

- Rows represent features.
- The first column contains feature identifiers.
- Every remaining column is a sample.
- Sample identifiers are the names of those remaining columns.
- Row names are not used as the feature-identifier input boundary.
- Values must represent already-preprocessed variance-stabilized expression.
- The package does not filter raw counts, normalize counts, or calculate VST
  values.
- Expression values are expected to be numeric.
- Existing behavior converts the expression matrix to numeric mode.
  Nonnumeric values may become `NA` with a warning.
- `NA`, `NaN`, and infinite values are not comprehensively rejected. They may
  propagate, be omitted by particular summaries, or be excluded from
  finite-score rankings according to the preserved calculation.
- Duplicate feature identifiers are not repaired or rejected.
- Duplicate sample-column identifiers are unsupported and must be avoided.
- Zero and negative VST values are not independently prohibited.
- Feature row order is preserved in the per-feature outputs and can determine
  the order of tied scores.
- At least one feature row and at least one named sample column are required.

## Metadata input

`metadata` must be a `data.frame` or `data.table` containing these columns:

- `sample_id`
- `condition`
- `baseline_block`
- `group_type`

The required values and normalization behavior are:

- `sample_id` is converted to character and surrounding whitespace is removed.
  It is not lowercased.
- `condition`, `baseline_block`, and `group_type` are converted to character,
  trimmed, and lowercased.
- After normalization and sample matching, `condition` must include the literal
  values `uv` and `des`.
- `group_type` uses the literal values `treatment` and `control`.
- The effect subtraction is treatment mean minus control mean.
- `baseline_block` identifies the subgroups used for stratification.
- Each condition must have at least two subgroup values.
- Subgroup labels may be arbitrary, but they are sorted lexicographically.
- If a condition has more than two subgroup values, only the first two
  lexicographically sorted values are used for `delta_1`, `delta_2`, subgroup
  weighting, and `delta_strat`.
- In that greater-than-two-subgroup case, `delta_pool` and `var_proxy` still use
  all matched samples in the condition. This preserved behavior can therefore
  make the pooled and stratified calculations use different subgroup scopes.
- Duplicate metadata sample identifiers are not repaired or rejected and must
  be avoided.
- Additional metadata columns are ignored.

### Treatment/control composition

Version 0.1.1 validates that matched samples include both `uv` and `des`, but it
does not validate treatment/control coverage.

To obtain complete subgroup and pooled effects, users should provide:

- at least one `treatment` and one `control` sample in each of the two
  lexicographically selected `baseline_block` subgroups; and
- consequently, at least one matched `treatment` and one matched `control`
  sample in each condition.

These composition requirements are not actively enforced. If a selected
subgroup lacks either role, its subgroup effect is `NA`, and the stratified
effect, Instability score, sign-flip result, and related summaries may also be
`NA` or `NaN`. If a condition lacks either role globally, `delta_pool` is
`NA`.

Subgroups beyond the first two lexicographically sorted values are technically
accepted. They are excluded from `delta_1`, `delta_2`, subgroup weighting, and
`delta_strat`, but all their matched samples remain part of `var_proxy`, and
their recognized treatment/control samples remain part of `delta_pool`.
Additional subgroups have no separately enforced treatment/control composition.
Users should avoid designs with more than two subgroups unless they explicitly
intend this preserved asymmetric behavior.

Samples with missing or unrecognized `group_type` values are excluded from
treatment and control means but remain included in `var_proxy`. Missing
expression values are not removed by the preserved `rowMeans()` and `var()`
calculations and can therefore produce missing feature-level results.

## Sample matching and ordering

Samples are matched using the intersection of `metadata$sample_id` and the VST
sample-column names.

- Samples present in only one input are excluded.
- Partial overlap is accepted.
- No overlap is an error.
- Both UV and DES must remain represented after matching.
- The VST and metadata do not need to start in the same order.
- Metadata is keyed by `sample_id` after matching.
- Inputs are not silently expanded to include unmatched samples.

Users should verify the matched design before analysis, particularly when
either input contains samples absent from the other.

## Fixed scientific behavior

The following choices are fixed and are not customizable through the public
API in version 0.1.1:

- conditions: `uv` and `des`;
- roles: `treatment` and `control`;
- treatment-minus-control subtraction;
- lexicographic subgroup sorting;
- first-two subgroup behavior;
- subgroup weights based on the treatment-plus-control sample count in each of
  the selected subgroups;
- sign epsilon: `1e-6`;
- Top-N: 150;
- tail K values: 50, 100, 200, 500, 1000, and 2000;
- headline K values: 200 and 1000;
- descending ranking by the relevant score;
- preserved source/current order for relevant ties;
- finite-score filtering in rank-displacement calculations;
- the existing Pearson and Spearman correlation behavior;
- the existing missing-value behavior of each calculation;
- private intermediate TSV serialization and rereading.

These settings are intentionally fixed to preserve the validated scientific
implementation. Exposing or changing them would require separate scientific
validation.

## How Instability is calculated

Within each condition:

- `delta_1` is the treatment mean minus control mean in subgroup 1.
- `delta_2` is the treatment mean minus control mean in subgroup 2.
- `delta_pool` is the treatment mean minus control mean across all matched
  samples in that condition.
- `delta_strat` is the subgroup-size-weighted average of `delta_1` and
  `delta_2`.
- `hetero_gap` is `abs(delta_1 - delta_2)`.
- `discordance` is `abs(delta_pool - delta_strat)`.
- `instability_score` is `hetero_gap + discordance`.

A `sign_flip` is `TRUE` when both subgroup effects exceed the fixed epsilon in
magnitude and have opposite signs.

Higher Instability scores identify features whose estimated effects are more
sensitive to subgroup structure or to pooling. They do not by themselves
establish causality, biological mechanism, or general validity outside the
locked UV/DES design.

## Runnable own-data example

This synthetic example uses arbitrary feature and sample identifiers. It does
not use frozen validation data.

```r
vst <- data.frame(
  feature = c("gene_A", "gene_B", "gene_C"),
  uv_a_control = c(1, 2, 3),
  uv_a_treatment = c(2, 1, 5),
  uv_b_control = c(4, 3, 2),
  uv_b_treatment = c(2, 6, 3),
  des_a_control = c(2, 3, 4),
  des_a_treatment = c(4, 2, 7),
  des_b_control = c(5, 4, 3),
  des_b_treatment = c(3, 7, 4),
  check.names = FALSE
)

metadata <- data.frame(
  sample_id = names(vst)[-1],
  condition = rep(c("uv", "des"), each = 4),
  baseline_block = rep(
    c("block_a", "block_a", "block_b", "block_b"),
    2
  ),
  group_type = rep(
    c("control", "treatment", "control", "treatment"),
    2
  )
)

result <- instability::compute_uv_des_instability(
  vst = vst,
  metadata = metadata
)

names(result)
result$uv_pool_vs_strat_all
result$des_pool_vs_strat_all
result$core_instability_signal_table
```

For an own dataset, replace the synthetic values and identifiers while
retaining the documented UV/DES, treatment/control, sample-matching, and
subgroup contract.

## Return value

The function returns an ordinary named list containing exactly eleven
`data.table` objects, in this order:

1. `uv_pool_vs_strat_all`
   One row per feature containing the complete UV feature-level calculation.

2. `uv_pool_vs_strat_top150`
   Up to 150 UV features ranked by descending `instability_score`.

3. `uv_pool_vs_strat_summary`
   One-row UV summary of feature count, sign flips, median scores, and the
   Instability-versus-variance correlation.

4. `des_pool_vs_strat_all`
   One row per feature containing the complete DES feature-level calculation.

5. `des_pool_vs_strat_top150`
   Up to 150 DES features ranked by descending `instability_score`.

6. `des_pool_vs_strat_summary`
   One-row DES summary.

7. `pool_vs_strat_metrics`
   One row per condition comparing pooled, stratified, and feature-divergence
   views.

8. `tail_jaccard_metrics`
   One row per condition and fixed K value containing overlaps among ranked
   tail sets and summaries within the top combined set.

9. `rank_displacement_pool_vs_strat`
   One row per finite-ranked feature and condition, comparing heterogeneity-gap
   and discordance ranks.

10. `tail_headline_summary`
    One row per condition containing selected K=200 and K=1000 tail metrics.

11. `core_instability_signal_table`
    Final one-row-per-condition combination of pool/strat metrics and selected
    tail summaries.

For datasets with fewer than 150 features, each Top-150 table contains all
available feature rows. For a tail K larger than the available feature count,
the available feature set is used.

## Important output columns

### Feature-level columns

- `feature`: identifier from the first VST column.
- `condition`: `uv` or `des`.
- `subgroup_1`, `subgroup_2`: selected lexicographically sorted subgroup labels.
- `delta_1`, `delta_2`: subgroup-specific treatment-minus-control effects.
- `delta_pool`: pooled treatment-minus-control effect using all matched samples
  in the condition.
- `delta_strat`: weighted average of `delta_1` and `delta_2`.
- `sign_flip`: whether the two non-negligible subgroup effects have opposite
  signs.
- `hetero_gap`: absolute difference between subgroup-specific effects.
- `discordance`: absolute difference between pooled and stratified effects.
- `instability_score`: `hetero_gap + discordance`.
- `var_proxy`: sample variance across all matched samples in the condition.

### Condition-summary columns

- `n_features`: number of feature rows.
- `n_sign_flip`: number of sign-flip features under the preserved missing-value
  rule.
- `frac_sign_flip` or `sign_flip_rate`: sign-flip proportion.
- `median_hetero_gap`, `median_discordance`, `median_instability`: medians of
  the named quantities.
- `cor_instability_vs_var`: Pearson correlation between `instability_score` and
  `var_proxy` using pairwise-complete observations.
- `effect_divergence`: `abs(delta_1 - delta_2)`; this is the same literal
  calculation as `hetero_gap` in the current implementation.
- `spearman_instability_vs_effect_divergence`: finite-pair Spearman
  correlation.
- `spearman_pool_vs_strat`: finite-pair Spearman correlation between pooled and
  stratified effects.
- `median_abs_pool_minus_strat`: median absolute pooled-minus-stratified
  difference.
- `q95_abs_pool_minus_strat`: 95th percentile of that absolute difference.

### Tail and rank columns

- `K`: fixed requested tail-set size.
- `j_gap_vs_dis`: Jaccard overlap between top `hetero_gap` and top
  `discordance` feature sets.
- `j_gap_vs_comb`: Jaccard overlap between top `hetero_gap` and top
  `instability_score` sets.
- `j_dis_vs_comb`: Jaccard overlap between top `discordance` and top
  `instability_score` sets.
- `flip_rate_in_top_comb`: sign-flip rate in the top Instability set.
- `median_abs_pool_minus_strat_in_top_comb`,
  `median_hetero_gap_in_top_comb`, and
  `median_discordance_in_top_comb`: medians within that set.
- `rank_strat`: descending rank by `hetero_gap`.
- `rank_pool`: descending rank by `discordance`.
- `signed_rank_diff`: `rank_pool - rank_strat`.
- `abs_rank_diff`: absolute rank difference.
- `x`: `log10(1 + abs_rank_diff)`.
- `j_gap_vs_dis_K200` and `j_gap_vs_dis_K1000`: selected Jaccard metrics.
- `flip_top200` and `flip_top1000`: selected top-set sign-flip rates.

Historically inherited column names are documented literally rather than
renamed.

## Scope and limitations

Version 0.1.1 does not provide:

- arbitrary condition selection;
- general subgroup handling;
- upstream raw-count filtering;
- normalization;
- VST generation;
- Gamma analysis;
- homogeneous-null controls;
- figure-generation workflows;
- annotation or NCBI workflows;
- synthetic benchmark workflows;
- archived or historical manuscript workflows.

The package does not repair duplicate identifiers, impute missing values, or
infer experimental roles. Inputs outside the documented contract may produce
missing results, warnings, or errors according to the preserved implementation.

Instability is a descriptive comparison of pooled and subgroup-aware effect
estimates. It is not a causal statistic and should not be interpreted as proof
of a biological mechanism.

## Provenance and citation

Authoritative scientific source repository:

https://github.com/sajadshahbaz/NM2_instability_signal

Authoritative scientific source commit:

`6a261c2aa665eb457a94aa1cead2cc8b39be1987`

Scientific-source baseline tag:

`v1.1.1`

Package repository:

https://github.com/sajadshahbaz/instability-r-package-dev

Scientific-source repository provenance DOI:

https://doi.org/10.5281/zenodo.21110170

That DOI identifies provenance for the scientific-source repository. It is not
a DOI for `instability` package version 0.1.1. No package DOI is claimed here.

The package citation metadata names Sajad Shahbazi as the software author.

## Validation

The packaged scientific implementation has previously reproduced all eleven
frozen UV/DES scientific reference tables byte-for-byte and passed exact parsed
and semantic comparison. Version 0.1.1 must repeat the complete equivalence
gate before release.

This validation does not extend to conditions or workflows outside the package
scope described above.
