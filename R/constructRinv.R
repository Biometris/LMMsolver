#' Helper function for constructing Rinv
#'
#' @keywords internal
constructRinv <- function(df,
                          residual,
                          weights) {
  ## If weights not defined, set equal to one.
  if (is.null(weights)) {
    weights <- rep(1, nrow(df))
  }
  if (!is.null(residual)) {
    lRinv <- list()
    column <- all.vars(residual)
    levels_f <- unique(df[[column]])
    cnt <- vector(length = length(levels_f))
    for (i in 1:length(levels_f)) {
      lRinv[[i]] <- spam::diag.spam(weights * (df[[column]] == levels_f[i]))
      cnt[i] <- sum(df[[column]] == levels_f[i])
    }
    names(lRinv) <- paste0(column, "_", levels_f, "!R")
    attr(lRinv, "cnt") <- cnt
  } else {
    lRinv <- list(residual = spam::diag.spam(weights))
    attr(lRinv, "cnt") <- nrow(df)
  }
  return(lRinv)
}
