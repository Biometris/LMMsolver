## some help functions, make part of LMM solve library:

#' @export
ndxMatrix <- function(df,
                      lZ,
                      Names) {
  n <- length(lZ)
  dim <- sapply(X = lZ, FUN = ncol)
  e <- cumsum(dim) + ncol(df)
  s <- e - dim + 1
  lM <- list()
  for (i in 1:n) {
    lM[[i]] <- c(s[i]:e[i])
  }
  names(lM) <- Names
  return(lM)
}

#' equally placed knots:
#'
#' @export
PsplinesKnots <- function(xmin,
                          xmax,
                          degree,
                          nseg) {
  dx <- (xmax - xmin) / nseg
  knots <- seq(xmin - degree * dx, xmax + degree * dx, by = dx)
  attr(knots, "degree") <- degree
  return(knots)
}

#' just to simplify splineDesign
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
