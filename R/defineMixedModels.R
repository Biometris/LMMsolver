# components of Rinv = phi1*A[[1]]+phi2*A[[2]]:
makeRlist <- function(df,
                      column) {
  lRinv <- list()
  levels_f <- levels(df[[column]])
  cnt <- vector(length = length(levels_f))
  for (i in 1:length(levels_f)) {
    lRinv[[i]] <- spam::diag.spam(1.0 * (df[[column]] == levels_f[i]))
    cnt[i] <- sum(df[[column]] == levels_f[i])
  }
  names(lRinv) <- paste0(column, "_", levels_f, "!R")
  attr(lRinv, "cnt") <- cnt

  return(lRinv)
}
