#' Construct Design matrix for predictions splines.
#'
#' @noRd
#' @keywords internal
designMatrixPredSplines <- function(object,
                                    grid = NULL,
                                    newdata = NULL,
                                    deriv = 0,
                                    includeIntercept = FALSE,
                                    which = 1) {
  ## Get dimension of fitted spline component.
  splRes <- object$splRes[[which]]
  splF_name <- splRes$term.labels.f
  splR_name <- splRes$term.labels.r
  ## Get content from splRes.
  x <- splRes$x
  knots <- splRes$knots
  scaleX <- splRes$scaleX
  pord <- splRes$pord
  degree <- splRes$degree
  splDim <- length(x)
  if (splDim == 1 && (!is.numeric(deriv) || length(deriv) > 1 || deriv < 0 ||
                      deriv != round(deriv))) {
    stop("deriv should be an integer greater than or equal to zero.\n")
  }
  if (splDim > 1 && deriv != 0) {
    deriv <- 0
    warning("deriv is ignored for ", splDim, "-dimensional splines.\n",
            call. = FALSE)
  }
  if (deriv > degree) {
    stop(deriv,
         "-order derivatives cannot be computed for B-splines of degree ",
         degree)
  }
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
    Bx <- mapply(FUN = Bsplines, knots, xGrid, deriv)
    ## Compute Bx over all dimensions.
    BxTot <- Reduce(RowKronecker, Bx)
  } else {
    if (!is.numeric(grid) || length(grid) != splDim) {
      stop("grid should be a numeric vector with length equal to the dimension ",
           "of the fitted spline: ", splDim,".\n")
    }
    ## Construct grid for each dimension.
    xGrid <- lapply(X = seq_len(splDim), FUN = function(i) {
      seq(attr(knots[[i]], which = 'xmin'), attr(knots[[i]], which = 'xmax'),
          length = grid[i])
    })
    ## Compute Bx per dimension.
    Bx <- mapply(FUN = Bsplines, knots, xGrid, deriv)
    ## Compute Bx over all dimensions.
    BxTot <- Reduce(`%x%`, Bx)
  }
  ## Compute G per dimension.
  G <- lapply(X=knots, FUN = function(x) {
    constructG(knots = x, scaleX = scaleX, pord = pord)})
  ## Compute G over all dimensions
  GTot <- Reduce('%x%', G)
  ## no scaling for first column of GTot
  GTot[,1] <- 1
  XTot <- BxTot %*% GTot
  ## Remove intercept (needed when fitting model to avoid singularities).
  XTot <- removeIntercept(XTot)

  labels <- c(object$term.labels.f, object$term.labels.r)
  lU <- list()
  dim <- object$dim
  for (i in seq_along(dim)) {
    lU[[i]] = spam::spam(x = 0, nrow = nrow(BxTot), ncol = dim[i])
  }
  if (includeIntercept) {
    lU[[1]] = spam::spam(x = 1, nrow = nrow(BxTot), ncol = 1)
  }
  if (!is.null(splRes$term.labels.f)) {
    ndx.f <- which(splRes$term.labels.f == labels)
    lU[[ndx.f]] <- XTot
  }
  ndx.r <- which(splRes$term.labels.r == labels)
  lU[[ndx.r]] <- BxTot

  U <- Reduce(spam::cbind.spam, lU)
  attr(U, "xGrid") <- xGrid
  return(U)
}

#' Obtain Smooth Trend.
#'
#' Obtain the smooth trend for models fitted with a spline component.
#'
#' @param object An object of class LMMsolve.
#' @param grid A numeric vector having the length of the dimension of the fitted
#' spline component. This represents the number of grid points at which a
#' surface will be computed.
#' @param newdata A data.frame containing new points for which the smooth
#' trend should be computed. Column names should include the names used when
#' fitting the spline model.
#' @param includeIntercept Should the value of the intercept be included in
#' the computed smooth trend? Ignored if deriv > 0.
#' @param deriv Derivative of B-splines, default 0. At the moment only
#' implemented for spl1D.
#' @param which An integer, for if there are multiple splxD terms in the model.
#' Default value is 1.
#'
#' @returns A data.frame with predictions for the smooth trend on the specified
#' grid. The standard errors are saved if `deriv` has default value 0.
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
#' ## Obtain the smooth trend for the fitted model on a dense grid.
#' smooth1 <- obtainSmoothTrend(LMM1_spline,
#'                             grid = 100)
#'
#' ## Obtain the smooth trend on a new data set - plots 10 to 40.
#' newdat <- data.frame(plot = 10:40)
#' smooth2 <- obtainSmoothTrend(LMM1_spline,
#'                             newdata = newdat)
#'
#' ## The first derivative of the smooth trend can be obtained by setting deriv = 1.
#' smooth3 <- obtainSmoothTrend(LMM1_spline,
#'                             grid = 100,
#'                             deriv = 1)
#'
#' ## For examples of higher order splines see the vignette.
#'
#' @export
obtainSmoothTrend <- function(object,
                              grid = NULL,
                              newdata = NULL,
                              deriv = 0,
                              includeIntercept = FALSE,
                              which = 1) {
  if (!inherits(object, "LMMsolve")) {
    stop("object should be an object of class LMMsolve.\n")
  }
  if (is.null(object$splRes)) {
    stop("The model was fitted without a spline component.\n")
  }
  if (is.null(grid) && is.null(newdata)) {
    stop("Specify either grid or newdata.\n")
  }
  if (!is.numeric(which) || which > length(object$splRes)) {
    stop("which should be an integer with value at most the number of fitted",
         "spline components.\n")
  }
  if (object$family$family == "multinomial") {
    stop("For multinomial response use predict function.\n")
  }

  splRes <- object$splRes[[which]]
  x <- splRes$x
  ## make the design matrix needed for predictions and corresponding
  ## standard errors.
  U <- designMatrixPredSplines(object, grid, newdata, deriv,
                               includeIntercept, which)
  xGrid <- attr(U,which="xGrid")
  ## calculate the predictions
  eta <- as.vector(U %*% object$coefMME)
  family <- object$family
  familyPred <- family$linkinv(eta)
  ## Construct output data.frame.
  if (!is.null(newdata)) {
    outDat <- newdata
    outDat[["ypred"]] <- familyPred
  } else {
    outDat <- data.frame(expand.grid(rev(xGrid)), ypred = familyPred)
    colnames(outDat)[-ncol(outDat)] <- rev(names(x))
    outDat <- outDat[c(names(x), "ypred")]
  }
  ## only add standard errors if:
  if (deriv == 0 && includeIntercept) {
    outDat[["se"]] <- calcStandardErrors(C = object$C, D = U)*abs(family$mu.eta(eta))
  }
  return(outDat)
}

