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

