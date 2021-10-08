#' fit 3D P-splines
#'
#' fit 3D P-splines using sparse implementation.
#'
#' @param x1 numerical vector containing the values of \code{x1} covariate.
#' @param x2 numerical vector containing the values of \code{x2} covariate.
#' @param x3 numerical vector containing the values of \code{x3} covariate.
#' @param nseg number of segments
#' @param scaleX logical, scale fixed effect or not. Default is TRUE,
#' no scaling.
#' @param pord order of penalty, default \code{pord=2}
#' @param degree degree of B-spline basis, default \code{degree=3}
#' @param x1lim numerical vector of length 2 containing the domain of covariate
#' \code{x1} where the knots should be placed. Default set to \code{NULL}
#' (covariate range).
#' @param x2lim numerical vector of length 2 containing the domain of covariate
#' \code{x2} where the knots should be placed. Default set to \code{NULL}
#' (covariate range).
#' @param x3lim numerical vector of length 2 containing the domain of covariate
#' \code{x3} where the knots should be placed. Default set to \code{NULL}
#' (covariate range).
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{X} - design matrix for fixed effect. The intercept is not included.
#'   \item \code{Z} - design matrix for random effect.
#'   \item \code{lGinv} - a list of precision matrices
#'   \item \code{knots} - a list of vectors with knot positions
#' }
#'
#' @importFrom stats setNames
#'
#' @export
spl3D <- function(x1,
                  x2,
                  x3,
                  nseg,
                  pord = 2,
                  degree = 3,
                  scaleX = TRUE,
                  x1lim = range(x1),
                  x2lim = range(x2),
                  x3lim = range(x3)) {
  ## Checks.
  if (!is.numeric(pord) || length(pord) > 1 || !pord %in% 1:2) {
    stop("pord should be either 1 or 2.\n")
  }
  if (!is.numeric(degree) || length(degree) > 1 || degree < 1 ||
      degree != round(degree)) {
    stop("degree should be a positive integer.\n")
  }
  if (!is.numeric(nseg) || length(nseg) != 3 || any(nseg < 1) ||
      any(nseg != round(nseg))) {
    stop("nseg should be a vector of length three of positive integers.\n")
  }
  ## Save names of the x-variables so they can be used later on in predictions.
  x1Name <- deparse(substitute(x1))
  x2Name <- deparse(substitute(x2))
  x3Name <- deparse(substitute(x3))
  xNames <- c(x1Name, x2Name, x3Name)
  missVars <- xNames[!sapply(X = xNames, FUN = exists,
                             where = parent.frame(), inherits = FALSE)]
  if (length(missVars) > 0) {
    stop("The following variables in the spline part of the model ",
         "are not in the data:\n", paste0(missVars, collapse = ", "), "\n",
         call. = FALSE)
  }
  if (!is.numeric(x1lim) || length(x1lim) != 2 ||
      x1lim[1] > range(x1)[1] || x1lim[2] < range(x1)[2]) {
    stop("x1lim should be a vector of length two with all values of ", x1Name,
         " between its lower and upper value.\n")
  }
  if (!is.numeric(x2lim) || length(x2lim) != 2 ||
      x2lim[1] > range(x2)[1] || x2lim[2] < range(x2)[2]) {
    stop("x2lim should be a vector of length two with all values of ", x2Name,
         " between its lower and upper value.\n")
  }
  if (!is.numeric(x3lim) || length(x3lim) != 2 ||
      x3lim[1] > range(x3)[1] || x3lim[2] < range(x3)[2]) {
    stop("x3lim should be a vector of length two with all values of ", x3Name,
         " between its lower and upper value.\n")
  }

  knots <- list(1)
  knots[[1]] <- PsplinesKnots(x1lim[1], x1lim[2], degree = degree, nseg = nseg[1])
  knots[[2]] <- PsplinesKnots(x2lim[1], x2lim[2], degree = degree, nseg = nseg[2])
  knots[[3]] <- PsplinesKnots(x3lim[1], x3lim[2], degree = degree, nseg = nseg[3])

  B1 <- Bsplines(knots[[1]], x1)
  B2 <- Bsplines(knots[[2]], x2)
  B3 <- Bsplines(knots[[3]], x3)

  q1 <- ncol(B1)
  q2 <- ncol(B2)
  q3 <- ncol(B3)

  # we have to calculate RowKronecker product only once:
  one.1 <- matrix(1, 1, ncol(B1))
  one.2 <- matrix(1, 1, ncol(B2))
  one.3 <- matrix(1, 1, ncol(B3))
  B123 <- (B1 %x% one.2 %x% one.3) *
          (one.1 %x% B2 %x% one.3) *
          (one.1 %x% one.2 %x% B3)

  DtD1 <- constructPenalty(q1, pord)
  DtD2 <- constructPenalty(q2, pord)
  DtD3 <- constructPenalty(q3, pord)

  X1 <- constructX(B1, x1, scaleX, pord)
  X2 <- constructX(B2, x2, scaleX, pord)
  X3 <- constructX(B3, x3, scaleX, pord)

  X <- RowKronecker(RowKronecker(X1, X2), X3)

  ## Remove intercept column to avoid singularity problems.
  X  <- removeIntercept(X)

  CCt1 <- constructCCt(q1, pord)
  CCt2 <- constructCCt(q2, pord)
  CCt3 <- constructCCt(q3, pord)

  CCt <- CCt1 %x% CCt2 %x% CCt3

  lGinv <- list()
  lGinv[[1]] <- DtD1 %x% spam::diag.spam(q2) %x% spam::diag.spam(q3) + CCt
  lGinv[[2]] <- spam::diag.spam(q1) %x% DtD2 %x% spam::diag.spam(q3) + CCt
  lGinv[[3]] <- spam::diag.spam(q1) %x% spam::diag.spam(q2) %x% DtD3 + CCt

  names(lGinv) <- c("s(x1)", "s(x2)", "s(x3)")

  if (is.null(X)) {
    dim.f <- NULL
    term.labels.f <- NULL
  } else {
    dim.f <- ncol(X)
    term.labels.f <- "splF"
  }
  dim.r <- ncol(B123)
  term.labels.r <- "splR"

  xList <- setNames(list(x1, x2, x3), c(x1Name, x2Name, x3Name))
  return(list(X = X, Z = B123, lGinv = lGinv, knots = knots,
              dim.f = dim.f, dim.r = dim.r, term.labels.f = term.labels.f,
              term.labels.r = term.labels.r, x = xList,
              pord = pord, degree = degree, scaleX = scaleX))
}
