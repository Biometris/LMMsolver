#' fit 1D P-splines
#'
#' fit 1D P-splines using sparse implementation.
#'
#' @param x numerical vector containing the values of x covariate.
#' @param nseg number of segments
#' @param scaleX logical, scale fixed effect or not. Default is FALSE,
#' no scaling.
#' @param pord order of penalty, default \code{pord=2}
#' @param degree degree of B-spline basis, default \code{degree=3}
#' @param xlim numerical vector of length 2 containing the domain of covariate
#' x where the knots should be placed. Default set to \code{NULL} (covariate range).
#'
#' @usage spl1D(x, nseg, pord=2, degree=3, scaleX=FALSE, xlim=NULL)
#'
#' @return A list with three matrices: fixed effect \code{X} (not including the
#' intercept), the random effect $Z$, and the precision matrix \code(Ginv}$ for
#' the random effect \code{Z}
#'
#' @export
spl1D <- function(x,
                  nseg,
                  pord = 2,
                  degree = 3,
                  scaleX = FALSE,
                  xlim = NULL) {
  if (is.null(xlim))
    xlim <- c(min(x), max(x))

  knots <- PsplinesKnots(xlim[1], xlim[2], degree = degree, nseg = nseg)

  B <- spam::as.spam(Bsplines(knots, x))
  q <- ncol(B)

  D <- spam::diff.spam(diag(q), diff = pord)
  DtD <- crossprod(D)

  if (pord == 2) {
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
    C <- spam::spam(x = 0, nrow = q, ncol = pord)
    C[1, 1] = C[q,2] = 1
    CCt <- spam::tcrossprod(C)
  } else {  # pord = 1
    X = NULL
    C <- spam::spam(x = 0, nrow = q, ncol = pord)
    C[1,1] = C[q,1] = 1
    # spam::tcrossprod doesn't work....
    CCt <- C %*% t(C)
  }

  lGinv <- list()
  lGinv[[1]] <- spam::as.spam(DtD + CCt)
  names(lGinv) <- "s(x)"

  return(list(X = X, Z = B, lGinv = lGinv))
}

