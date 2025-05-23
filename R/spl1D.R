#' Fit P-splines
#'
#' Fit multi dimensional P-splines using sparse implementation.
#'
#' @param x,x1,x2,x3 The variables in the data containing the values of
#' the \code{x} covariates.
#' @param nseg The number of segments
#' @param scaleX Should the fixed effects be scaled.
#' @param pord The order of penalty, default \code{pord = 2}
#' @param degree The degree of B-spline basis, default \code{degree = 3}
#' @param cyclic Cyclic or linear B-splines; default \code{cyclic=FALSE}
#' @param xlim,x1lim,x2lim,x3lim A numerical vector of length 2 containing the
#' domain of the corresponding x covariate where the knots should be placed.
#' Default set to \code{NULL}, when the covariate range will be used.
#' @param cond Conditional factor: splines are defined conditional on the level.
#' Default \code{NULL}.
#' @param level The level of the conditional factor. Default \code{NULL}.
#'
#' @returns A list with the following elements:
#' \itemize{
#'   \item \code{X} - design matrix for fixed effect. The intercept is not included.
#'   \item \code{Z} - design matrix for random effect.
#'   \item \code{lGinv} - a list of precision matrices
#'   \item \code{knots} - a list of vectors with knot positions
#'   \item \code{dim.f} - the dimensions of the fixed effect.
#'   \item \code{dim.r} - the dimensions of the random effect.
#'   \item \code{term.labels.f} - the labels for the fixed effect terms.
#'   \item \code{term.labels.r} - the labels for the random effect terms.
#'   \item \code{x} - a list of vectors for the spline variables.
#'   \item \code{pord} - the order of the penalty.
#'   \item \code{degree} - the degree of the B-spline basis.
#'   \item \code{scaleX} - logical indicating if the fixed effects are scaled.
#'   \item \code{EDnom} - the nominal effective dimensions.
#' }
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit a model with a 1-dimensional spline at the plot level.
#' LMM1_spline <- LMMsolve(fixed = yield ~ rep + gen,
#'                        spline = ~spl1D(x = plot, nseg = 20),
#'                        data = john.alpha)
#'
#' summary(LMM1_spline)
#'
#' ## Fit model on US precipitation data from spam package.
#' data(USprecip, package = "spam")
#'
#' ## Only use observed data
#' USprecip <- as.data.frame(USprecip)
#' USprecip <- USprecip[USprecip$infill == 1, ]
#'
#' ## Fit a model with a 2-dimensional P-spline.
#' LMM2_spline <- LMMsolve(fixed = anomaly ~ 1,
#'                        spline = ~spl2D(x1 = lon, x2 = lat, nseg = c(41, 41)),
#'                        data = USprecip)
#'
#' summary(LMM2_spline)
#'
#' @seealso \code{\link{LMMsolve}}
#'
#' @importFrom stats setNames
#'
#' @export
spl1D <- function(x,
                  nseg,
                  pord = 2,
                  degree = 3,
                  cyclic = FALSE,
                  scaleX = TRUE,
                  xlim = range(x),
                  cond = NULL,
                  level = NULL) {
  ## Checks
  if (!is.numeric(pord) || length(pord) > 1 || !pord %in% 1:3) {
    stop("pord should be equal to 1, 2 or 3.\n")
  }
  if (!is.numeric(degree) || length(degree) > 1 || degree < 1 ||
      degree != round(degree)) {
    stop("degree should be a positive integer.\n")
  }
  if (pord > degree) {
    stop("pord should be less or equal to degree.\n")
  }
  if (!is.numeric(nseg) || length(nseg) > 1 || nseg < 1 ||
      nseg != round(nseg)) {
    stop("nseg should be a positive integer.\n")
  }
  if (!is.logical(cyclic) || length(cyclic) > 1) {
    stop("cyclic should be FALSE or TRUE.\n")
  }
  if (cyclic) {
    xlim = c(0,1)
    if (min(x) < 0 || max(x) > 1) {
      stop("x should be in the range [0,1] for cyclic data.\n")
    }
    if (pord != 2) {
      stop("pord should be equal to two for cyclic data.\n")
    }
  }

  ## Save names of the x-variables so they can be used later on in predictions.
  xName <- deparse(substitute(x))
  if (!exists(xName, where = parent.frame(), inherits = FALSE)) {
    stop("The following variables in the spline part of the model ",
         "are not in the data:\n", xName, "\n",
         call. = FALSE)
  }
  ## check (syntax) conditional factor.
  conditional <- checkConditionalFactor(x, cond, level)
  if (conditional) {
    Nelem <- length(x)
    ndx <- cond == level
    x <- x[ndx]
  }
  checkLim(lim = xlim, limName = "xlim", x = x, xName = xName)
  knots <- vector(mode = "list", length = 1)
  knots[[1]] <- PsplinesKnots(xlim[1], xlim[2], degree = degree, nseg = nseg, cyclic=cyclic)
  B <- Bsplines(knots[[1]], x)
  q <- ncol(B)
  if (conditional) {
    sel <- which(ndx)
    B <- extSpamMatrix(B, sel, length(ndx))
  }
  G <- constructG(knots[[1]], scaleX, pord)
  X <- B %*% G
  ## nominal effective dimension.
  EDnom <- ncol(B) - ncol(X)
  ## Remove intercept column to avoid singularity problems.
  X <- removeIntercept(X)
  ## Construct list of sparse precision matrices.
  scaleFactor <- calcScaleFactor(knots, pord)
  lGinv <- constructGinvSplines(q, knots, pord, scaleFactor)
  names(lGinv) <- paste0("s(", xName, ")")
  if (is.null(X)) {
    dim.f <- NULL
    term.labels.f <- NULL
  } else {
    dim.f <- ncol(X)
    if (pord > 2) {
      term.labels.f <- paste0("pol(", xName, ")")
    } else {
      term.labels.f <- paste0("lin(", xName, ")")
    }
  }
  dim.r <- ncol(B)
  term.labels.r <- paste0("s(", xName, ")")
  if (conditional) {
    if (!is.null(term.labels.f)) {
      term.labels.f <- paste0(term.labels.f, "_", level)
    }
    term.labels.r <- paste0(term.labels.r, "_", level)
    names(lGinv) <- paste0("s(", xName, ")_", level)
  }
  xList <- setNames(list(x), xName)
  return(list(X = X, Z = B, lGinv = lGinv, knots = knots,
              dim.f = dim.f, dim.r = dim.r, term.labels.f = term.labels.f,
              term.labels.r = term.labels.r, x = xList, pord = pord,
              degree = degree, scaleX = scaleX, EDnom = EDnom,
              scaleFactor = scaleFactor))
}

