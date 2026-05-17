#' Create a data frame for use in making predictions.
#'
#' Constructs a grid of values spanning the ranges of the spline covariates
#' stored in an \code{LMMsolver} object.
#'
#' @param object An \code{LMMsolver} object.
#' @param grid A numeric vector specifying the number of grid points for each
#'   spline dimension. Its length must equal the number of spline variables.
#'
#' @return A data frame containing all combinations of grid values for the
#'   spline covariates.
#'
#' @details
#' For each spline variable, equally spaced values are generated between the
#' minimum and maximum values of the B-splines. The Cartesian product of these sequences
#' is returned using \code{\link[base]{expand.grid}}.
#'
#' @examples
#' \dontrun{
#' ## Create a 200 x 300 grid for a two-dimensional spline term
#' grd <- makeGrid(fit, grid = c(200, 300))
#'
#' head(grd)
#' }
#'
#' @export
makeGrid <- function(object, grid) {
  splRes <- object$splRes
  if (is.null(splRes)) {
    stop("Spline not defined.\n")
  }
  if (length(splRes) > 1) {
    stop("Not implemented yet: multiple spline terms.\n")
  }
  if (!is.numeric(grid)) {
    stop("grid should be a numeric vector.\n")
  }
  spl <- splRes[[1]]
  if (length(grid) != length(spl$knots)) {
    stop("Argument dim has the wrong length.\n")
  }
  VarNames <- names(spl$x)
  knots <- spl$knots
  xmin_val <- sapply(knots, FUN = function(x) {
    attr(x, which = "xmin")
  })
  xmax_val <- sapply(knots, FUN = function(x) {
    attr(x, which = "xmax")
  })
  dim <- length(knots)

  L <- list()
  for (i in seq_len(dim)) {
    L[[i]] <- seq(xmin_val[i], xmax_val[i], length.out = grid[i])
  }
  df <- expand.grid(L)
  colnames(df) <- VarNames
  return(df)
}


