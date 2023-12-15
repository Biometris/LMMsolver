#' @describeIn spl1D 3-dimensional splines
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
  checkLim(lim = x1lim, limName = "x1lim", x = x1, xName = x1Name)
  checkLim(lim = x2lim, limName = "x2lim", x = x2, xName = x2Name)
  checkLim(lim = x3lim, limName = "x3lim", x = x3, xName = x3Name)
  knots <- vector(mode = "list", length = 3)
  knots[[1]] <- PsplinesKnots(x1lim[1], x1lim[2], degree = degree, nseg = nseg[1])
  knots[[2]] <- PsplinesKnots(x2lim[1], x2lim[2], degree = degree, nseg = nseg[2])
  knots[[3]] <- PsplinesKnots(x3lim[1], x3lim[2], degree = degree, nseg = nseg[3])
  B1 <- Bsplines(knots[[1]], x1)
  B2 <- Bsplines(knots[[2]], x2)
  B3 <- Bsplines(knots[[3]], x3)
  q <- c(ncol(B1), ncol(B2), ncol(B3))
  B123 <- RowKronecker(RowKronecker(B1, B2), B3)
  G1 <- constructG(knots[[1]], scaleX, pord)
  G2 <- constructG(knots[[2]], scaleX, pord)
  G3 <- constructG(knots[[3]], scaleX, pord)
  G <- G1 %x% G2 %x% G3
  X <- B123 %*% G
  ## nominal effective dimension.
  EDnom = rep(ncol(B123) - ncol(X), 3)
  ## Remove intercept column to avoid singularity problems.
  X <- removeIntercept(X)
  ## Construct list of sparse precision matrices.
  scaleFactor <- calcScaleFactor(knots, pord)
  lGinv <- constructGinvSplines(q, knots, pord, scaleFactor)
  names(lGinv) <- paste0("s(", xNames, ")")
  if (is.null(X)) {
    dim.f <- NULL
    term.labels.f <- NULL
  } else {
    dim.f <- ncol(X)
    term.labels.f <- paste0("lin(", paste(xNames, collapse = ", "), ")")
  }
  dim.r <- ncol(B123)
  term.labels.r <- paste0("s(", paste(xNames, collapse = ", "), ")")
  xList <- setNames(list(x1, x2, x3), xNames)
  return(list(X = X, Z = B123, lGinv = lGinv, knots = knots,
              dim.f = dim.f, dim.r = dim.r, term.labels.f = term.labels.f,
              term.labels.r = term.labels.r, x = xList,
              pord = pord, degree = degree, scaleX = scaleX, EDnom = EDnom,
              scaleFactor = scaleFactor))
}
