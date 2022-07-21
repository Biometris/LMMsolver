#' @describeIn spl1D 2-dimensional splines
#'
#' @export
spl2D <- function(x1,
                  x2,
                  nseg,
                  pord = 2,
                  degree = 3,
                  scaleX = TRUE,
                  x1lim = range(x1),
                  x2lim = range(x2)) {
  ## Checks.
  if (!is.numeric(pord) || length(pord) > 1 || !pord %in% 1:2) {
    stop("pord should be either 1 or 2.\n")
  }
  if (!is.numeric(degree) || length(degree) > 1 || degree < 1 ||
      degree != round(degree)) {
    stop("degree should be a positive integer.\n")
  }
  if (!is.numeric(nseg) || length(nseg) != 2 || any(nseg < 1) ||
      any(nseg != round(nseg))) {
    stop("nseg should be a vector of length two of positive integers.\n")
  }
  ## Save names of the x-variables so they can be used later on in predictions.
  x1Name <- deparse(substitute(x1))
  x2Name <- deparse(substitute(x2))
  xNames <- c(x1Name, x2Name)
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
  knots <- list()
  knots[[1]] <- PsplinesKnots(x1lim[1], x1lim[2], degree = degree, nseg = nseg[1])
  knots[[2]] <- PsplinesKnots(x2lim[1], x2lim[2], degree = degree, nseg = nseg[2])
  B1 <- Bsplines(knots[[1]], x1)
  B2 <- Bsplines(knots[[2]], x2)
  q <- c(ncol(B1), ncol(B2))
  B12 <- RowKronecker(B1, B2)
  X1 <- constructX(B1, x1, scaleX, pord)
  X2 <- constructX(B2, x2, scaleX, pord)
  X <- RowKronecker(X1, X2)
  ## nominal effective dimension.
  EDnom = rep(ncol(B12) - ncol(X), 2)
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
  dim.r <- ncol(B12)
  term.labels.r <-  paste0("s(", paste(xNames, collapse = ", "), ")")
  xList <- setNames(list(x1, x2), xNames)
  return(list(X = X, Z = B12, lGinv = lGinv, knots = knots,
              dim.f = dim.f, dim.r = dim.r, term.labels.f = term.labels.f,
              term.labels.r = term.labels.r, x = xList, pord = pord,
              degree = degree, scaleX = scaleX, EDnom = EDnom, scaleFactor=scaleFactor))

}

