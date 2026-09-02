.prepare_uv_des_inputs <- function(vst, metadata) {
  if (!inherits(vst, "data.frame")) {
    stop("vst must be a data.frame or data.table", call. = FALSE)
  }
  if (!inherits(metadata, "data.frame")) {
    stop("metadata must be a data.frame or data.table", call. = FALSE)
  }

  required <- c("sample_id", "condition", "baseline_block", "group_type")
  missing <- setdiff(required, names(metadata))
  if (length(missing) > 0) {
    stop(
      sprintf("metadata missing required columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  meta <- data.table::copy(data.table::as.data.table(metadata))
  vst_dt <- data.table::copy(data.table::as.data.table(vst))

  sid_col <- "sample_id"
  cond_col <- "condition"
  sub_col <- "baseline_block"
  grp_col <- "group_type"

  meta[, (sid_col) := trimws(as.character(get(sid_col)))]
  meta[, (cond_col) := tolower(trimws(as.character(get(cond_col))))]
  meta[, (sub_col) := tolower(trimws(as.character(get(sub_col))))]
  meta[, (grp_col) := tolower(trimws(as.character(get(grp_col))))]

  features <- vst_dt[[1]]
  vst_samples <- names(vst_dt)[-1]

  common <- intersect(meta[[sid_col]], vst_samples)
  if (length(common) == 0) {
    stop("No sample overlap between meta and VST.", call. = FALSE)
  }
  meta <- meta[get(sid_col) %in% common]
  data.table::setkeyv(meta, sid_col)

  X <- as.matrix(vst_dt[, ..common])
  mode(X) <- "numeric"
  rownames(X) <- features

  list(X = X, meta = meta)
}

.effect_delta <- function(mat, treat_ids, ctrl_ids) {
  mt <- if (length(treat_ids) > 0) rowMeans(mat[, treat_ids, drop = FALSE]) else rep(NA_real_, nrow(mat))
  mc <- if (length(ctrl_ids) > 0) rowMeans(mat[, ctrl_ids, drop = FALSE]) else rep(NA_real_, nrow(mat))
  mt - mc
}

.sign_with_epsilon <- function(x, eps = 1e-6) {
  if (is.na(x)) return(NA_integer_)
  if (abs(x) <= eps) return(0L)
  if (x > 0) return(1L)
  -1L
}

.compute_uv_des_condition <- function(X, meta, condition) {
  cond <- condition
  sid_col <- "sample_id"
  cond_col <- "condition"
  sub_col <- "baseline_block"
  grp_col <- "group_type"

  m <- meta[get(cond_col) == cond]
  if (nrow(m) == 0) {
    return(invisible(NULL))
  }

  subs <- sort(unique(m[[sub_col]]))
  if (length(subs) < 2) {
    stop(sprintf("Condition '%s' has <2 subgroups in %s.", cond, sub_col), call. = FALSE)
  }
  if (length(subs) > 2) {
    subs <- subs[1:2]
  }
  s1 <- subs[1]
  s2 <- subs[2]

  treat_all <- m[get(grp_col) == "treatment", get(sid_col)]
  ctrl_all <- m[get(grp_col) == "control", get(sid_col)]

  m1 <- m[get(sub_col) == s1]
  m2 <- m[get(sub_col) == s2]

  treat_1 <- m1[get(grp_col) == "treatment", get(sid_col)]
  ctrl_1 <- m1[get(grp_col) == "control", get(sid_col)]
  treat_2 <- m2[get(grp_col) == "treatment", get(sid_col)]
  ctrl_2 <- m2[get(grp_col) == "control", get(sid_col)]

  d1 <- .effect_delta(X, treat_1, ctrl_1)
  d2 <- .effect_delta(X, treat_2, ctrl_2)
  d_pool <- .effect_delta(X, treat_all, ctrl_all)

  n1 <- length(treat_1) + length(ctrl_1)
  n2 <- length(treat_2) + length(ctrl_2)
  w1 <- n1 / (n1 + n2)
  w2 <- n2 / (n1 + n2)
  d_strat <- w1 * d1 + w2 * d2

  cond_samples <- m[[sid_col]]
  var_all <- apply(X[, cond_samples, drop = FALSE], 1, var)

  eps <- 1e-6
  sgn1 <- vapply(d1, .sign_with_epsilon, integer(1), eps = eps)
  sgn2 <- vapply(d2, .sign_with_epsilon, integer(1), eps = eps)

  sign_flip <- (sgn1 != 0L) & (sgn2 != 0L) & (sgn1 != sgn2)
  hetero_gap <- abs(d1 - d2)
  discordance <- abs(d_pool - d_strat)
  instability_score <- hetero_gap + discordance

  data.table::data.table(
    feature = rownames(X),
    condition = cond,
    subgroup_1 = s1,
    subgroup_2 = s2,
    delta_1 = d1,
    delta_2 = d2,
    delta_pool = d_pool,
    delta_strat = d_strat,
    sign_flip = sign_flip,
    hetero_gap = hetero_gap,
    discordance = discordance,
    instability_score = instability_score,
    var_proxy = var_all
  )
}

.make_top150 <- function(condition_table) {
  out_rank <- condition_table[order(-instability_score)]
  out_rank[1:min(150L, .N)]
}

.summarize_uv_des_condition <- function(condition_table) {
  condition_table[, .(
    n_features = .N,
    n_sign_flip = sum(sign_flip, na.rm = TRUE),
    frac_sign_flip = mean(sign_flip, na.rm = TRUE),
    median_hetero_gap = median(hetero_gap, na.rm = TRUE),
    median_discordance = median(discordance, na.rm = TRUE),
    median_instability = median(instability_score, na.rm = TRUE),
    cor_instability_vs_var = suppressWarnings(cor(instability_score, var_proxy, use = "pairwise.complete.obs"))
  )]
}
