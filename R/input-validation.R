.validate_uv_des_structure <- function(vst, metadata) {
  if (inherits(vst, "data.frame")) {
    if (ncol(vst) == 0L) {
      stop(
        "vst must contain a first feature-identifier column and at least one sample column",
        call. = FALSE
      )
    }
    if (nrow(vst) == 0L) {
      stop("vst must contain at least one feature row", call. = FALSE)
    }
    if (ncol(vst) < 2L) {
      stop("vst must contain at least one named sample column", call. = FALSE)
    }
  }

  if (inherits(metadata, "data.frame") && nrow(metadata) == 0L) {
    stop("metadata must contain at least one sample row", call. = FALSE)
  }

  invisible(NULL)
}

.validate_required_uv_des_conditions <- function(prepared) {
  observed <- unique(prepared$meta$condition)
  missing <- setdiff(c("uv", "des"), observed)
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "metadata must contain matched samples for required condition(s): %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}
