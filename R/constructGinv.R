constructGinv <- function(ginverse, lGinv, dim.r, term.labels.r) {
  ## Replace identity matrix with ginverse for random terms.
  if (!is.null(ginverse)) {
    ## Convert to spam.
    e <- cumsum(dim.r)
    s <- e - dim.r + 1
    for (i in seq_along(ginverse)) {
      ginvVar <- names(ginverse)[i]
      ginvMat <- ginverse[[i]]
      k <- which(term.labels.r == ginvVar)
      if (length(k) == 0) {
        stop("ginverse element ", ginvVar, " not defined in random part.\n", call. = FALSE)
      }
      if (dim.r[k] != nrow(ginvMat)) {
        stop("Dimensions of ", ginvVar, " should match number of levels ",
             "for corresponding factor in data.\n", call. = FALSE)
      }
      ndx <- s[k]:e[k]
      ## as.spam drops row and column names, so only converting to spam
      ## after all checking is done.
      lGinv[[k]][ndx, ndx] <- spam::as.spam(ginvMat)
    }
  }
  return(lGinv)
}
