#' sap2D, without spectral decomposition
#'
#' sap2D, without spectral decomposition
#'
#' @param x ...
#' @param nseg ...
#' @param scaleX ...
#' @param pord ...
#' @param degr ...
#'
#' @export
spl1D <- function(x,
                  nseg,
                  pord = 2,
                  degr = 3,
                  scaleX = FALSE) {
  xlim <- c(min(x) - 0.01, max(x) + 0.01)

  knots <- PsplinesKnots(xlim[1], xlim[2], degree = degr, nseg = nseg)

  B <- spam::as.spam(Bsplines(knots, x))
  q <- ncol(B)

  D <- spam::diff.spam(diag(q), diff = pord)
  DtD <- crossprod(D)

  if (scaleX) {
    ## calculate the linear/fixed parts.
    U_null <- cbind(1, scale(1:q))

    U_null <- apply(U_null, MARGIN = 2, function(x) (x / normVec(x)))

    X <- B %*% U_null
    ## Remove intercept column to avoid singularity problems.
    X <- X[, -1, drop=FALSE]
  } else {
    X <- cbind(1, x)
    ## Remove intercept column to avoid singularity problems.
    X <- X[, -1, drop=FALSE]
  }
  # assuming pord=2
  C <- spam::spam(x = 0, nrow = q, ncol = pord)
  C[1, 1] = C[q,2] = 1

  CCt <- spam::tcrossprod(C)

  lGinv <- list()
  lGinv[[1]] <- spam::as.spam(DtD + CCt)
  names(lGinv) <- c("x")

  return(list(X = X, Z = B, lGinv = lGinv))
}

