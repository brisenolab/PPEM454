#' Summarize pairwise identity in an alignment
#' @param alignment Character vector of equal-length aligned sequences.
#' @return List with pairwise percentage identity matrix and overall, min,
#'   and max percentages. Overall is the unweighted mean across unique pairs.
#' @details Each pair uses only columns where both symbols are A, C, G, T,
#'   or U (U is treated as T). Gaps and ambiguity symbols are excluded.
#'   Pairs with no comparable sites are NA and excluded from summaries.
#' @export
sequence_similarity <- function(alignment) {
  if (!is.character(alignment) || length(alignment) < 2L || anyNA(alignment) ||
      any(!nzchar(alignment)) || length(unique(nchar(alignment))) != 1L) {
    stop("Supply at least two nonempty, equal-length aligned sequences.")
  }
  chars <- strsplit(chartr("U", "T", toupper(alignment)), "", fixed = TRUE)
  ids <- names(alignment)
  if (is.null(ids)) ids <- paste0("sequence_", seq_along(chars))
  result <- matrix(NA_real_, length(chars), length(chars), dimnames = list(ids, ids))
  for (i in seq_along(chars)) for (j in seq_len(i)) {
    valid <- chars[[i]] %in% c("A", "C", "G", "T") & chars[[j]] %in% c("A", "C", "G", "T")
    value <- if (any(valid)) 100 * mean(chars[[i]][valid] == chars[[j]][valid]) else NA_real_
    result[i, j] <- result[j, i] <- value
  }
  pairs <- result[upper.tri(result)]
  pairs <- pairs[!is.na(pairs)]
  list(pairwise = result, overall = if (length(pairs)) mean(pairs) else NA_real_,
       min = if (length(pairs)) min(pairs) else NA_real_,
       max = if (length(pairs)) max(pairs) else NA_real_)
}
