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

# some help functions for spl1D, spl2D and spl3D functions:

constructPenalty <-function(q, pord)
{
  D <- spam::diff.spam(diag(q), diff = pord)
  DtD <- crossprod(D)
}

constructX <- function(B, x, scaleX, pord)
{
  if (pord == 2) {
    if (scaleX) {
      ## calculate the linear/fixed parts.
      U_null <- cbind(1, scale(1:q))
      U_null <- apply(U_null, MARGIN = 2, function(x) (x / normVec(x)))
      X <- B %*% U_null
    } else {
      X <- cbind(1, x)
    }
  } else {  # pord = 1
    X = matrix(data=1, nrow=length(x), ncol=1)
  }
  X
}

constructConstraint <- function(q, pord)
{
  if (pord == 2) {
    C <- spam::spam(x = 0, nrow = q, ncol = pord)
    C[1, 1] = C[q,2] = 1
  } else {  # pord = 1
    C <- spam::spam(x = 0, nrow = q, ncol = pord)
    C[1,1] = C[q,1] = 1
    # spam::tcrossprod doesn't work....
  }
  C
}

removeIntercept <- function(X) {
  if (ncol(X)==1) {
    X = NULL
  } else {
    X <- X[, -1,  drop=FALSE]
  }
  X
}


