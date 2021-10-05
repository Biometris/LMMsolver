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
