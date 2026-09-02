\name{compute_uv_des_instability}
\alias{compute_uv_des_instability}
\title{Compute the UV and Desiccation Instability Outputs}
\usage{
compute_uv_des_instability(vst, metadata)
}
\arguments{
  \item{vst}{A \code{data.frame} or \code{data.table}. The first column
  contains feature identifiers and every remaining column is a named sample.
  Values must represent already-preprocessed variance-stabilized expression.}

  \item{metadata}{A \code{data.frame} or \code{data.table} containing the
  canonical columns \code{sample_id}, \code{condition},
  \code{baseline_block}, and \code{group_type}.}
}
\value{
An ordinary named list containing exactly eleven \code{data.table} objects,
in this order:
\itemize{
  \item \code{uv_pool_vs_strat_all}
  \item \code{uv_pool_vs_strat_top150}
  \item \code{uv_pool_vs_strat_summary}
  \item \code{des_pool_vs_strat_all}
  \item \code{des_pool_vs_strat_top150}
  \item \code{des_pool_vs_strat_summary}
  \item \code{pool_vs_strat_metrics}
  \item \code{tail_jaccard_metrics}
  \item \code{rank_displacement_pool_vs_strat}
  \item \code{tail_headline_summary}
  \item \code{core_instability_signal_table}
}
}
\description{
Computes the validated Instability outputs for the existing ultraviolet (UV)
and desiccation (DES) conditions from an already-preprocessed VST expression
table and its associated metadata.
}
\details{
The function is limited to the existing UV/DES implementation. It does not
perform raw-count filtering, normalization, variance-stabilizing
transformation, Gamma analysis, null-control analysis, annotation, or
plotting.

Scientifically relevant intermediate table serialization and rereading occur
inside a private temporary directory. The temporary directory is removed when
the function exits. No output-directory argument is exposed, and the function
does not write scientific outputs into the caller working directory.
}
