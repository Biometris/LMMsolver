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
  if (!is.numeric(pord) || length(pord) > 1 || !pord %in% 1:2) {
    stop("pord should be either 1 or 2.\n")
  }
  if (!is.numeric(degree) || length(degree) > 1 || degree < 1 ||
      degree != round(degree)) {
    stop("degree should be a positive integer.\n")
  }

  xName <- deparse(substitute(x))

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
  X <- removeIntercept(X)

  lGinv <- list()
  lGinv[[1]] <- spam::as.spam(DtD + CCt)
  names(lGinv) <- "s(x)"

  if (is.null(X)) {
    dim.f <- NULL
    term.labels.f <- NULL
  } else {
    dim.f <- c(ncol(X))
    term.labels.f <- "splF"
  }
  dim.r <- ncol(B)
  term.labels.r <- "splR"

  xList <- setNames(list(x), xName)
  return(list(X = X, Z = B, lGinv = lGinv, knots = knots,
              dim.f = dim.f, dim.r = dim.r, term.labels.f = term.labels.f,
              term.labels.r = term.labels.r, x = xList, pord = pord,
              degree = degree, scaleX = scaleX))
}




