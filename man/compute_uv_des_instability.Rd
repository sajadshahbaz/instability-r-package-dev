\name{compute_uv_des_instability}
\alias{compute_uv_des_instability}
\title{Compute the UV and Desiccation Instability Outputs}
\usage{
compute_uv_des_instability(vst, metadata)
}
\arguments{
  \item{vst}{A \code{data.frame} or \code{data.table} containing
  already-preprocessed variance-stabilized expression data. Rows are features.
  The first column contains feature identifiers, and every remaining named
  column is a sample whose column name is its sample identifier. Row names are
  not used as the feature-identifier input boundary. Values are expected to be
  numeric. Nonnumeric values may be converted to \code{NA} by the preserved
  numeric coercion. Missing and nonfinite values are not comprehensively
  rejected and may propagate or be excluded according to each calculation.
  Duplicate feature identifiers are not repaired; duplicate sample-column
  identifiers are unsupported and must be avoided. At least one feature row
  and one named sample column are required.}

  \item{metadata}{A \code{data.frame} or \code{data.table} containing
  \code{sample_id}, \code{condition}, \code{baseline_block}, and
  \code{group_type}. Sample identifiers are converted to character and
  trimmed. Condition, subgroup, and group-role values are converted to
  character, trimmed, and lowercased. After matching, \code{condition} must
  include literal \code{uv} and \code{des}; \code{group_type} uses literal
  \code{treatment} and \code{control}. Each condition requires at least two
  \code{baseline_block} values. More than two subgroup values are accepted
  under the preserved first-two behavior described in Details.}
}
\value{
An ordinary named list containing exactly eleven \code{data.table} objects:
\describe{
  \item{\code{uv_pool_vs_strat_all}}{One row per feature containing the
  complete UV feature-level calculation.}
  \item{\code{uv_pool_vs_strat_top150}}{Up to 150 UV features ranked by
  descending Instability score.}
  \item{\code{uv_pool_vs_strat_summary}}{One-row UV aggregate summary.}
  \item{\code{des_pool_vs_strat_all}}{One row per feature containing the
  complete DES feature-level calculation.}
  \item{\code{des_pool_vs_strat_top150}}{Up to 150 DES features ranked by
  descending Instability score.}
  \item{\code{des_pool_vs_strat_summary}}{One-row DES aggregate summary.}
  \item{\code{pool_vs_strat_metrics}}{One row per condition comparing pooled,
  stratified, and feature-divergence views.}
  \item{\code{tail_jaccard_metrics}}{One row per condition and fixed K value
  containing ranked-tail overlaps and top-set summaries.}
  \item{\code{rank_displacement_pool_vs_strat}}{One row per finite-ranked
  feature and condition comparing heterogeneity-gap and discordance ranks.}
  \item{\code{tail_headline_summary}}{One row per condition containing selected
  K=200 and K=1000 tail metrics.}
  \item{\code{core_instability_signal_table}}{The final one-row-per-condition
  combination of pool/strat metrics and selected tail summaries.}
}
}
\description{
Computes the validated Instability outputs for the existing ultraviolet (UV)
and desiccation (DES) conditions from an already-preprocessed VST expression
table and associated metadata.
}
\details{
The function implements the locked UV/DES design. It does not support arbitrary
condition names.

Samples are matched using the intersection of \code{metadata$sample_id} and the
VST sample-column names. Samples present in only one input are excluded.
Partial overlap is accepted; no overlap is an error. Both UV and DES must remain
represented after matching. Metadata is keyed by sample identifier after
matching.

Within each condition, subgroup labels are sorted lexicographically. At least
two subgroup values are required. If more than two exist, only the first two
are used for subgroup-specific effects, subgroup weighting, and the stratified
effect. The pooled effect and variance proxy continue to use all matched
samples in that condition.

Version 0.1.1 validates the presence of matched \code{uv} and \code{des}
samples, but it does not validate treatment/control coverage. To obtain
complete subgroup and pooled effects, each of the two selected subgroups should
contain at least one \code{treatment} and one \code{control} sample. This also
ensures that each condition contains both roles.

These composition requirements are not actively enforced. If a selected
subgroup lacks either role, its subgroup effect is \code{NA}, and the
stratified effect, Instability score, sign-flip result, and related summaries
may be \code{NA} or \code{NaN}. If a condition lacks either role globally,
\code{delta_pool} is \code{NA}.

Additional subgroups beyond the first two are technically accepted. They are
excluded from \code{delta_1}, \code{delta_2}, subgroup weighting, and
\code{delta_strat}. Their recognized treatment/control samples remain part of
\code{delta_pool}, and all their matched samples remain part of
\code{var_proxy}. No treatment/control composition is separately enforced for
these additional subgroups. Users should avoid more than two subgroups unless
they explicitly intend this preserved asymmetric behavior.

Samples with missing or unrecognized \code{group_type} values are excluded
from treatment and control means but remain included in \code{var_proxy}.
Missing expression values are not removed by the preserved \code{rowMeans()}
and \code{var()} calculations and can produce missing feature-level results.

The effect columns have these literal meanings:
\itemize{
  \item \code{delta_1}: subgroup 1 treatment mean minus control mean.
  \item \code{delta_2}: subgroup 2 treatment mean minus control mean.
  \item \code{delta_pool}: treatment mean minus control mean across all matched
  samples in the condition.
  \item \code{delta_strat}: subgroup-size-weighted average of
  \code{delta_1} and \code{delta_2}; subgroup sizes count treatment plus
  control samples.
  \item \code{hetero_gap}: \code{abs(delta_1 - delta_2)}.
  \item \code{discordance}: \code{abs(delta_pool - delta_strat)}.
  \item \code{instability_score}: \code{hetero_gap + discordance}.
  \item \code{sign_flip}: both subgroup effects exceed the fixed epsilon in
  magnitude and have opposite signs.
  \item \code{var_proxy}: sample variance across all matched samples in the
  condition.
}

Important summary and ranking columns include:
\itemize{
  \item \code{rank_strat}: descending rank by \code{hetero_gap}.
  \item \code{rank_pool}: descending rank by \code{discordance}.
  \item \code{signed_rank_diff}: \code{rank_pool - rank_strat}.
  \item \code{abs_rank_diff}: absolute rank difference.
  \item \code{x}: \code{log10(1 + abs_rank_diff)}.
  \item \code{K}: fixed requested tail-set size.
  \item \code{j_gap_vs_dis}: Jaccard overlap between top
  \code{hetero_gap} and top \code{discordance} sets.
  \item \code{j_gap_vs_comb}: Jaccard overlap between top
  \code{hetero_gap} and top \code{instability_score} sets.
  \item \code{j_dis_vs_comb}: Jaccard overlap between top
  \code{discordance} and top \code{instability_score} sets.
  \item \code{flip_rate_in_top_comb}: sign-flip rate within the top
  Instability set.
  \item \code{effect_divergence}: \code{abs(delta_1 - delta_2)}, the same
  literal calculation as \code{hetero_gap} in the current implementation.
  \item \code{cor_instability_vs_var}: Pearson correlation using
  pairwise-complete observations.
  \item \code{spearman_instability_vs_effect_divergence} and
  \code{spearman_pool_vs_strat}: finite-pair Spearman correlations.
}

The following choices are fixed and are not public parameters: conditions
\code{uv} and \code{des}; roles \code{treatment} and \code{control};
treatment-minus-control subtraction; lexicographic first-two subgroup
selection; subgroup-count weighting; sign epsilon \code{1e-6}; Top-N
\code{150}; tail K values \code{50}, \code{100}, \code{200}, \code{500},
\code{1000}, and \code{2000}; headline K values \code{200} and \code{1000};
descending ranking; preserved tie behavior; finite-score rank filtering;
correlation methods; missing-value behavior; and private TSV serialization and
rereading.

For K values larger than the available feature count, all available features
are used. Relevant score ties retain the preserved source/current ordering.

The package does not perform raw-count filtering, normalization, VST
generation, Gamma analysis, null-control analysis, figure generation,
annotation, or synthetic benchmark workflows. It does not repair duplicate
identifiers, impute missing values, or infer experimental roles.

Higher Instability scores identify features whose effects are more sensitive
to subgroup structure or pooling. Instability is descriptive; it is not a
causal statistic or proof of a biological mechanism.
}
\examples{
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

result <- compute_uv_des_instability(vst, metadata)

names(result)
result$uv_pool_vs_strat_all
result$des_pool_vs_strat_all
result$core_instability_signal_table
}
