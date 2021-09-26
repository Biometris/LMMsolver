#' fit 1D P-splines
#'
#' fit 1D P-splines using sparse implementation.
#'
#' @param x numerical vector containing the values of x covariate.
#' @param nseg number of segments
#' @param scaleX logical, scale fixed effect or not. Default is TRUE,
#' no scaling.
#' @param pord order of penalty, default \code{pord=2}
#' @param degree degree of B-spline basis, default \code{degree=3}
#' @param xlim numerical vector of length 2 containing the domain of covariate
#' x where the knots should be placed. Default set to \code{NULL} (covariate range).
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{X} - design matrix for fixed effect. The intercept is not included.
#'   \item \code{Z} - design matrix for random effect.
#'   \item \code{lGinv} - a list of precision matrices
#'   \item \code{knots} - a list of vectors with knot positions
#' }
#'
#' @export
spl1D <- function(x,
                  nseg,
                  pord = 2,
                  degree = 3,
                  scaleX = TRUE,
                  xlim = NULL) {
  if (is.null(xlim))
    xlim <- c(min(x), max(x))

  knots <- list()
  knots[[1]] <- PsplinesKnots(xlim[1], xlim[2], degree = degree, nseg = nseg)

  B <- spam::as.spam(Bsplines(knots[[1]], x))
  q <- ncol(B)

  DtD <- constructPenalty(q, pord)

  X <- constructX(B, x, scaleX, pord)
  CCt <- constructCCt(q, pord)

  ## Remove intercept column to avoid singularity problems.
  X  <- removeIntercept(X)

  lGinv <- list()
  lGinv[[1]] <- spam::as.spam(DtD + CCt)
  names(lGinv) <- "s(x)"

  if (is.null(X))
  {
    dim.f = NULL
    term.labels.f = NULL
  } else {
    dim.f = c(ncol(X))
    term.labels.f = c('splF')
  }
  dim.r = c(ncol(B))
  term.labels.r = c('splR')
  return(list(X = X, Z = B, lGinv = lGinv, knots = knots,
              dim.f=dim.f, dim.r=dim.r, term.labels.f=term.labels.f,
              term.labels.r=term.labels.r, x=x, pord=pord, degree=degree,
              scaleX=scaleX))
}

#' obtain Smooth Trend for 1D P-splines
#'
#' @param obj an object of class LMMsolve
#' @param grid number of grid points.
#' @export
obtainSmoothTrend1D <- function(object, grid) {
  x <- object$splRes$x
  knots <- object$splRes$knots[[1]]

  xgrid <- seq(min(x), max(x), length = grid)

  Bx <- spam::as.spam(Bsplines(knots, xgrid))

  X <- constructX(Bx, xgrid, object$splRes$scaleX,object$splRes$pord)
  X <- removeIntercept(X)

  mu <- coef.LMMsolve(object)$'(Intercept)'
  if (is.null(X))
  {
    bc <- 0.0
  } else
  {
    bc <- as.vector(X %*% coef.LMMsolve(object)$splF)
  }
  sc <- as.vector(Bx %*% coef.LMMsolve(object)$splR)
  fit <- mu + bc + sc
  p.data <- list(x=xgrid)
  L <- list(p.data = p.data, eta=fit, mu=fit)
  return(L)
}


