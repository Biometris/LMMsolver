#' Row-wise kronecker product
#'
#' Row-wise kronecker product
#'
#' @param X1 A matrix.
#' @param X2 A matrix.
#'
#' @return The row-wise kronecker product of X1 and X2.
#'
#' @keywords internal
RowKronecker <- function(X1,
                         X2) {
  one.1 <- matrix(1, 1, ncol(X1))
  one.2 <- matrix(1, 1, ncol(X2))
  rowKron <- kronecker(X1, one.2) * kronecker(one.1, X2)
  return(rowKron)
}

#' Construct index matrix
#'
#' Construct index matrix.
#'
#' @param df A data.frame.
#' @param lZ A list of matrices.
#' @param names A character vector of names.
#'
#' @return The index matrix.
#'
#' @keywords internal
ndxMatrix <- function(df,
                      lZ,
                      names) {
  n <- length(lZ)
  dim <- sapply(X = lZ, FUN = ncol)
  e <- cumsum(dim) + ncol(df)
  s <- e - dim + 1
  lM <- list()
  for (i in 1:n) {
    lM[[i]] <- c(s[i]:e[i])
  }
  names(lM) <- names
  return(lM)
}

#' Construct equally placed knots
#'
#' Construct equally placed knots.
#'
#' @param xmin A numerical value.
#' @param xmax A numerical value.
#' @param degree A numerical value.
#' @param nseg A numerical value.
#'
#' @return A numerical vector of knot positions.
#'
#' @keywords internal
PsplinesKnots <- function(xmin,
                          xmax,
                          degree,
                          nseg) {
  dx <- (xmax - xmin) / nseg
  knots <- seq(xmin - degree * dx, xmax + degree * dx, by = dx)
  attr(knots, "degree") <- degree
  attr(knots, "xmin") <- xmin
  attr(knots, "xmax") <- xmax

  return(knots)
}

#' Construct design matrix for B-Splines
#'
#' Construct design matrix for B-Splines.
#'
#' @param knots A numerical vector of knot positions.
#' @param x a numeric vector of values at which to evaluate the B-spline
#' functions or derivatives.
#' @param deriv A numerical value. The derivative of the given order is
#' evaluated at the x positions.
#'
#' @export
Bsplines <- function(knots,
                     x,
                     deriv = 0) {
  degree <- attr(knots, "degree")
  B <- splines::splineDesign(knots = knots, x = x, ord = degree + 1,
                             derivs = deriv)
  return(B)
}
