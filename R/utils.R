#' Norm of a vector
#'
#' @param x a numerical vector
#' @return norm of vector x
#'
normVec <- function(x) {
  return(sqrt(sum(x ^ 2)))
}

#' combine two lists of Ginv matrices
#'
#' @param lGinv1 a list of sparse matrices
#' @param lGinv2 a list of sparse matrices
#' @description
#'
#' output is a combined list, of dimension (d1+d2) x (d1+d2),
#' with matrices in lGinv1 extended to [A1 0]
#'                                     [0  0]
#' and              lGinv2 extended to [0  0]
#'                                     [0 A2]
#' @return a list of sparse matrices of dimension (d1+d2) x (d1+d2)
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

#' construct P-splines penalty matrix D'D
#'
#' @param q integer with dimensions
#' @param pord order of the penalty
#'
#' @return qxq matrix D'D of class spam
#' @export
constructPenalty <-function(q, pord)
{
  D <- spam::diff.spam(diag(q), diff = pord)
  DtD <- crossprod(D)
  DtD
}

#' construct fixed part of the spline model
#'
#' @param B matrix with B-spline basis.
#' @param x vector with values for x
#' @param scaleX logical. If scaleX is FALSE, the original x is used. If scaleX is TRUE,
#' scaling is used, based on the B-splines basis. For details see the code
#' @param pord order of the penalty, values 1 or 2.
#'
#' @return a matrix X
#' @export
constructX <- function(B, x, scaleX, pord)
{
  q <- ncol(B)
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

#' construct constraint matrix
#'
#' @param q dimension of the B-spline basis used
#' @param pord order of the penalty matrix (pord=1 or 2).
#'
#' @return a q x q matrix of type spam
#' @export
constructCCt <- function(q, pord)
{
  if (pord == 2) {
    C <- spam::spam(x = 0, nrow = q, ncol = pord)
    C[1, 1] = C[q,2] = 1
  } else {  # pord = 1
    C <- spam::spam(x = 0, nrow = q, ncol = pord)
    C[1,1] = C[q,1] = 1
  }
  CCt <- C %*% t(C)
  CCt
}

#' remove the intercept from a design matrix
#'
#' @param X design matrix
#' @return a matrix if \code{X} has more than one column, otherwise return NULL
#' @export
removeIntercept <- function(X) {
  if (ncol(X)==1) {
    X = NULL
  } else {
    X <- X[, -1,  drop=FALSE]
  }
  X
}


