#' Row-wise kronecker product
#'
#' Row-wise kronecker product
#'
#' @param X1 A matrix.
#' @param X2 A matrix.
#'
#' @returns The row-wise kronecker product of X1 and X2.
#'
#' @keywords internal
RowKronecker <- function(X1,
                         X2) {
  if (inherits(X1, "spam") && inherits(X2, "spam")) {
    rowKron <- RowKron(X1, X2)
  } else {
    one.1 <- matrix(1, 1, ncol(X1))
    one.2 <- matrix(1, 1, ncol(X2))
    rowKron <- kronecker(X1, one.2) * kronecker(one.1, X2)
  }
  return(rowKron)
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
#' @returns A numerical vector of knot positions.
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
  attr(knots, "dx") <- dx

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
#' @keywords internal
Bsplines <- function(knots,
                     x,
                     deriv = 0) {
  degree <- attr(knots, "degree")
  B <- splines::splineDesign(knots = knots, x = x, ord = degree + 1,
                             derivs = deriv, sparse = TRUE, outer.ok = TRUE)
  return(spam::as.spam.dgCMatrix(B))
}
