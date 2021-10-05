#' Norm of a vector
#'
#' @param x a numerical vector.
#'
#' @return norm of vector x
#'
#' @keywords internal
normVec <- function(x) {
  return(sqrt(sum(x ^ 2)))
}

#' Combine two lists of Ginv matrices
#'
#' output is a combined list, of dimension (d1+d2) x (d1+d2),
#' with matrices in lGinv1 extended to [A1 0]
#'                                     [0  0]
#' and              lGinv2 extended to [0  0]
#'                                     [0 A2]
#' @return a list of sparse matrices of dimension (d1+d2) x (d1+d2)
#'
#' @param lGinv1 a list of sparse matrices.
#' @param lGinv2 a list of sparse matrices.
#'
#' @keywords internal
expandGinv <- function(lGinv1, lGinv2) {
  if (is.null(lGinv1)) return(lGinv2)
  if (is.null(lGinv2)) return(lGinv1)

  d1 <- nrow(lGinv1[[1]])
  d2 <- nrow(lGinv2[[1]])

  zero1 <- spam::diag.spam(0, d2)
  zero2 <- spam::diag.spam(0, d1)

  lGinv1a <- lapply(X = lGinv1, FUN = function(x) {
    spam::bdiag.spam(x, zero1)
  })
  lGinv2a <- lapply(X = lGinv2, FUN = function(x) {
    spam::bdiag.spam(zero2, x)
  })
  lGinv <- c(lGinv1a, lGinv2a)
  return(lGinv)
}

#' Construct P-splines penalty matrix D'D
#'
#' @param q integer with dimensions.
#' @param pord order of the penalty.
#'
#' @return qxq matrix D'D of class spam
#'
#' @keywords internal
constructPenalty <-function(q,
                            pord) {
  D <- spam::diff.spam(spam::diag.spam(q), diff = pord)
  DtD <- spam::crossprod.spam(D)
  return(DtD)
}

#' Construct fixed part of the spline model
#'
#' @param B matrix with B-spline basis.
#' @param x vector with values for x.
#' @param scaleX logical. If scaleX is FALSE, the original x is used. If scaleX
#' is TRUE, scaling is used, based on the B-splines basis. For details see
#' the code.
#' @param pord order of the penalty, values 1 or 2.
#'
#' @return a matrix X
#'
#' @keywords internal
constructX <- function(B,
                       x,
                       scaleX,
                       pord) {
  if (pord == 2) {
    if (scaleX) {
      ## calculate the linear/fixed parts.
      U_null <- cbind(1, scale(seq_len(ncol(B))))
      U_null <- apply(U_null, MARGIN = 2, function(x) (x / normVec(x)))
      X <- B %*% U_null
    } else {
      X <- cbind(1, x)
    }
  } else {  # pord = 1
    X <- matrix(data = 1, nrow = length(x), ncol = 1)
  }
  X
}

#' Construct constraint matrix
#'
#' @param q dimension of the B-spline basis used.
#' @param pord order of the penalty matrix (pord=1 or 2).
#'
#' @return a q x q matrix of type spam
#'
#' @keywords internal
constructCCt <- function(q,
                         pord) {
  C <- spam::spam(x = 0, nrow = q, ncol = pord)
  C[1, 1] <- C[q, pord] <- 1
  CCt <- C %*% t(C)
  return(CCt)
}

#' Remove the intercept from a design matrix
#'
#' @param X design matrix.
#'
#' @return a matrix if \code{X} has more than one column, otherwise return NULL
#'
#' @keywords internal
removeIntercept <- function(X) {
  if (ncol(X) == 1) {
    X <- NULL
  } else {
    X <- X[, -1,  drop = FALSE]
  }
  return(X)
}

#' Calculate the Nominal Effective dimension
#'
#' @keywords internal
calcNomEffDim <- function(X,
                          Z,
                          dim.r) {
  if (is.null(Z)) return(NULL)
  p <- ncol(X)

  EDnom <- vector(length = length(dim.r))
  e <- cumsum(dim.r)
  s <- e - dim.r + 1
  # for each variance component in Z:
  for (i in 1:length(dim.r)) {
    if (dim.r[i] > 100) {
      EDnom[i] <- dim.r[i]
    } else {
      ndx <- s[i]:e[i]
      XZ <- cbind(X, Z[, ndx])
      r <- qr(XZ)$rank
      EDnom[i] <- r - p
    }
  }
  return(EDnom)
}

#' Check variables in formula
#'
#' Check that all variables in a formula are present in the data.
#'
#' @importFrom utils hasName
#'
#' @keywords internal
checkFormVars <- function(formula,
                          data) {
  if (!is.null(formula)) {
    ## Get variables in formula.
    formVars <- all.vars(formula)
    missVars <- formVars[!hasName(data, formVars)]
    if (length(missVars) > 0) {
      ## Get the name of formula as passed to the function.
      formName <- deparse(substitute(formula))
      stop("The following variables in the ", formName, " part of the model ",
           "are not in the data:\n", paste0(missVars, collapse = ", "), "\n",
           call. = FALSE)
    }
  }
}

