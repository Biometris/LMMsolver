#' Obtain Smooth Trend.
#'
#' Obtain the smooth trend for models fitted with a spline component.
#'
#' @param object An object of class LMMsolve.
#' @param grid A numeric vector having the length of the dimension of the fitted
#' spline component. This represents the number of grid points at which a
#' surface will be computed.
#' @param newdata ...
#' @param includeIntercept Should the value of the intercept be included in
#' the computed smooth trend?
#'
#' @return A data.frame with predictions for the smooth trend on the specified
#' grid.
#'
#' @export
obtainSmoothTrend <- function(object,
                              grid = NULL,
                              newdata = NULL,
                              includeIntercept = FALSE) {
  if (!inherits(object, "LMMsolve")) {
    stop("object should be an object of class LMMsolve.\n")
  }
  if (is.null(object$splRes)) {
    stop("The model was fitted without a spline component.\n")
  }
  if (is.null(grid) && is.null(newdata)) {
    stop("Specify either grid or newdata.\n")
  }
  ## Get dimension of fitted spline component.
  splRes <- object$splRes
  ## Get content from splRes.
  x <- splRes$x
  knots <- splRes$knots
  scaleX <- splRes$scaleX
  pord <- splRes$pord
  if (!is.null(newdata)) {
    if (!inherits(newdata, "data.frame")) {
      stop("newdata should be a data.frame.\n")
    }
    missX <- names(x)[!sapply(X = names(x), FUN = function(name) {
      hasName(x = newdata, name = name)
    })]
    if (length(missX) > 0) {
      stop("The following smoothing variables are not in newdata:\n",
           paste0(missX, collapse = ", "), "\n")
    }
    ## Construct grid for each dimension.
    xGrid <- lapply(X = seq_along(x), FUN = function(i) {
      newdata[[names(x)[i]]]
    })
    ## Compute Bx per dimension.
    Bx <- mapply(FUN = Bsplines, knots, xGrid)
    ## Compute Bx over all dimensions.
    BxTot <- Reduce(RowKronecker, Bx)
  } else {
    splDim <- length(x)
    if (!is.numeric(grid) || length(grid) != splDim) {
      stop("grid should be a numeric vector with length equal to the dimension ",
           "of the fitted spline: ", splDim,".\n")
    }
    ## Construct grid for each dimension.
    xGrid <- lapply(X = seq_along(x), FUN = function(i) {
      seq(min(x[[i]]), max(x[[i]]), length = grid[i])
    })
    ## Compute Bx per dimension.
    Bx <- mapply(FUN = Bsplines, knots, xGrid)
    ## Compute Bx over all dimensions.
    BxTot <- Reduce(`%x%`, Bx)
  }
  ## Compute X per dimension.
  X <- mapply(FUN = function(x, y) {
    constructX(B = x, x = y, scaleX = scaleX, pord = pord)
  }, Bx, xGrid, SIMPLIFY = FALSE)
  ## Compute X over all dimensions.
  if (!is.null(newdata)) {
    XTot <- Reduce(RowKronecker, X)
  } else {
    XTot <- Reduce(`%x%`, X)
  }
  ## Remove intercept (needed when fitting model to avoid singularities).
  XTot <- removeIntercept(XTot)
  ## Get intercept and compute contribution of fixed and random terms.
  if (includeIntercept) {
    mu <- coef.LMMsolve(object)$'(Intercept)'
  } else {
    mu <- 0
  }
  if (is.null(XTot)) {
    bc <- 0
  } else {
    bc <- as.vector(XTot %*% coef.LMMsolve(object)$splF)
  }
  sc <- as.vector(BxTot %*% coef.LMMsolve(object)$splR)
  ## Compute fitted values.
  fit <- mu + bc + sc
  ## Construct output data.frame.
  if (!is.null(newdata)) {
    outDat <- newdata
    outDat[["ypred"]] <- fit
  } else {
    outDat <- data.frame(expand.grid(rev(xGrid)), ypred = fit)
    colnames(outDat)[-ncol(outDat)] <- rev(names(x))
    outDat <- outDat[c(names(x), "ypred")]
  }
  return(outDat)
}






#### Original functions below for comparison



#' Obtain Smooth Trend for 1D P-splines
#'
#' @param object an object of class LMMsolve.
#' @param grid number of grid points.
#'
#' @keywords internal
obtainSmoothTrend1D <- function(object,
                                grid) {
  x <- object$splRes$x[[1]]
  knots <- object$splRes$knots[[1]]

  xgrid <- seq(min(x), max(x), length = grid)

  Bx <- spam::as.spam(Bsplines(knots, xgrid))

  X <- constructX(Bx, xgrid, object$splRes$scaleX,object$splRes$pord)
  X <- removeIntercept(X)

  mu <- coef.LMMsolve(object)$'(Intercept)'
  if (is.null(X)) {
    bc <- 0
  } else {
    bc <- as.vector(X %*% coef.LMMsolve(object)$splF)
  }
  sc <- as.vector(Bx %*% coef.LMMsolve(object)$splR)
  fit <- mu + bc + sc
  p.data <- list(x = xgrid)
  L <- list(p.data = p.data, eta = fit, mu = mu)
  return(L)
}

#' Obtain Smooth Trend for 2D P-splines.
#'
#' @param object an object of class LMMsolve
#' @param grid a numeric vector of length 2, with the number of grid points at
#' which a two-dimensional surface will be computed.
#'
#' @keywords internal
obtainSmoothTrend2D <- function(object,
                                grid) {
  x1 <- object$splRes$x[[1]]
  x2 <- object$splRes$x[[2]]
  knots1 <- object$splRes$knots[[1]]
  knots2 <- object$splRes$knots[[2]]

  x1grid <- seq(min(x1), max(x1), length = grid[1])
  x2grid <- seq(min(x2), max(x2), length = grid[2])

  Bx1 <- spam::as.spam(Bsplines(knots1, x1grid))
  Bx2 <- spam::as.spam(Bsplines(knots2, x2grid))

  B12x <- Bx1 %x% Bx2

  X1 <- constructX(Bx1, x1grid, object$splRes$scaleX,object$splRes$pord)
  X2 <- constructX(Bx2, x2grid, object$splRes$scaleX,object$splRes$pord)
  X <- X1 %x% X2
  X <- removeIntercept(X)

  mu <- coef.LMMsolve(object)$'(Intercept)'
  if (is.null(X)) {
    bc <- 0
  } else {
    bc <- as.vector(X %*% coef.LMMsolve(object)$splF)
  }
  sc <- as.vector(B12x %*% coef.LMMsolve(object)$splR)
  fit <- mu + bc + sc

  p.data <- list(x1=x1grid, x2=x2grid)
  L <- list(p.data = p.data, eta = fit, mu = mu)

  return(L)
}

#' Obtain Smooth Trend for 3D P-splines
#'
#' @param object an object of class LMMsolve
#' @param grid a numeric vector of length 3, with the number of grid points at
#' which a three-dimensional surface will be computed.
#'
#' @keywords internal
obtainSmoothTrend3D <- function(object,
                                grid) {
  x1 <- object$splRes$x[[1]]
  x2 <- object$splRes$x[[2]]
  x3 <- object$splRes$x[[3]]

  knots1 <- object$splRes$knots[[1]]
  knots2 <- object$splRes$knots[[2]]
  knots3 <- object$splRes$knots[[3]]

  x1grid <- seq(min(x1), max(x1), length = grid[1])
  x2grid <- seq(min(x2), max(x2), length = grid[2])
  x3grid <- seq(min(x3), max(x3), length = grid[3])

  Bx1 <- spam::as.spam(Bsplines(knots1, x1grid))
  Bx2 <- spam::as.spam(Bsplines(knots2, x2grid))
  Bx3 <- spam::as.spam(Bsplines(knots3, x3grid))

  B123x <- Bx1 %x% Bx2 %x% Bx3

  X1 <- constructX(Bx1, x1grid, object$splRes$scaleX, object$splRes$pord)
  X2 <- constructX(Bx2, x2grid, object$splRes$scaleX, object$splRes$pord)
  X3 <- constructX(Bx3, x3grid, object$splRes$scaleX, object$splRes$pord)

  X <- X1 %x% X2 %x% X3
  X <- removeIntercept(X)

  mu <- coef.LMMsolve(object)$'(Intercept)'
  if (is.null(X)) {
    bc <- 0
  } else {
    bc <- as.vector(X %*% coef.LMMsolve(object)$splF)
  }
  sc <- as.vector(B123x %*% coef.LMMsolve(object)$splR)
  fit <- mu + bc + sc

  p.data <- list(x1 = x1grid, x2 = x2grid, x3 = x3grid)
  L <- list(p.data = p.data, eta = fit, mu = mu)

  return(L)
}
