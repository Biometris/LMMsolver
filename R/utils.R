#' @keywords internal
normVec <- function(x) {
  return(sqrt(sum(x ^ 2)))
}

#
# combine two lists of lGinv matrices
# lGinv1 are sparse matrices with dimension d1 x d1
# lGinv2 are sparse matrices with dimension d2 x d2
#
# output is a combined list, of dimension (d1+d2) x (d1+d2),
# with matrices in lGinv1 extended to [A1 0]
#                                     [0  0]
#
# and              lGinv2 extended to [0  0]
#                                     [0 A2]
ExpandGinv <- function(lGinv1, lGinv2)
{
  if (is.null(lGinv1)) return(lGinv2)
  if (is.null(lGinv2)) return(lGinv1)

  l1 <- length(lGinv1)
  l2 <- length(lGinv2)
  d1 <- nrow(lGinv1[[1]])
  d2 <- nrow(lGinv2[[1]])

  zero1 <- spam::diag.spam(0.0, d2)
  zero2 <- spam::diag.spam(0.0, d1)

  lGinv <- list()
  for (i in 1:l1) {
    lGinv[[i]] = spam::bdiag.spam(lGinv1[[i]], zero1)
  }
  for (i in 1:l2) {
    lGinv[[l1+i]] = spam::bdiag.spam(zero2, lGinv2[[i]])
  }
  names(lGinv) <- c(names(lGinv1), names(lGinv2))
  lGinv
}
