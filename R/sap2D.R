#' fit 2D P-splines
#'
#' fit 2D P-splines using sparse implementation.
#'
#' @param x1 numerical vector containing the values of \code{x1} covariate.
#' @param x2 numerical vector containing the values of \code{x2} covariate.

#' @param nseg number of segments
#' @param scaleX logical, scale fixed effect or not. Default is TRUE,
#' no scaling.
#' @param pord order of penalty, default \code{pord=2}
#' @param degree degree of B-spline basis, default \code{degree=3}
#' @param x1lim numerical vector of length 2 containing the domain of covariate
#' \code{x1} where the knots should be placed. Default set to \code{NULL} (covariate range).
#' @param x2lim numerical vector of length 2 containing the domain of covariate
#' \code{x2} where the knots should be placed. Default set to \code{NULL} (covariate range).
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{X} - design matrix for fixed effect. The intercept is not included.
#'   \item \code{Z} - design matrix for random effect.
#'   \item \code{lGinv} - a list of precision matrices
#'   \item \code{knots} - a list of vectors with knot positions
#' }
#'
#' @export
sap2D <- function(x1,
                  x2,
                  nseg,
                  pord = 2,
                  degree = 3,
                  scaleX = TRUE,
                  x1lim = NULL,
                  x2lim = NULL) {
  if (is.null(x1lim))
    x1lim <- c(min(x1), max(x1))
  if (is.null(x2lim))
    x2lim <- c(min(x2), max(x2))

  knots <- list()
  knots[[1]] <- PsplinesKnots(x1lim[1], x1lim[2], degree = degree, nseg = nseg[1])
  knots[[2]] <- PsplinesKnots(x2lim[1], x2lim[2], degree = degree, nseg = nseg[2])

  B1 <- spam::as.spam(Bsplines(knots[[1]], x1))
  B2 <- spam::as.spam(Bsplines(knots[[2]], x2))
  q1 <- ncol(B1)
  q2 <- ncol(B2)

  one.1 <- matrix(1, 1, ncol(B1))
  one.2 <- matrix(1, 1, ncol(B2))
  B12 <- (B1 %x% one.2) * (one.1 %x% B2)

  DtD1 <- constructPenalty(q1, pord)
  DtD2 <- constructPenalty(q2, pord)

  X1 <- constructX(B1, x1, scaleX, pord)
  X2 <- constructX(B2, x2, scaleX, pord)
  X <- RowKronecker(X1, X2)

  ## Remove intercept column to avoid singularity problems.
  X  <- removeIntercept(X)

  CCt1 <- constructCCt(q1, pord)
  CCt2 <- constructCCt(q2, pord)
  CCt <- CCt1 %x% CCt2


  lGinv <- list()
  lGinv[[1]] <- DtD1 %x% spam::diag.spam(q2) + CCt
  lGinv[[2]] <- spam::diag.spam(q1) %x% DtD2 + CCt

  names(lGinv) <- c("s(x1)", "s(x2)")

  if (is.null(X))
  {
    dim.f = NULL
    term.labels.f = NULL
  } else {
    dim.f = c(ncol(X))
    term.labels.f = c('splF')
  }
  dim.r = c(ncol(B12))
  term.labels.r = c('splR')

  return(list(X = X, Z = B12, lGinv = lGinv, knots = knots,
              dim.f=dim.f, dim.r=dim.r, term.labels.f=term.labels.f,
              term.labels.r=term.labels.r, x1=x1, x2=x2, pord=pord, degree=degree,
              scaleX=scaleX))

}

#' obtain Smooth Trend for 2D P-splines
#'
#' @param obj an object of class LMMsolve
#' @param grid a numeric vector of length 2, with the number of grid points at which a
#' two-dimensional surface will be computed.
#' @export
obtainSmoothTrend2D <- function(object, grid) {
  x1 <- object$splRes$x1
  x2 <- object$splRes$x2
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

  mu <- coef(object)$'(Intercept)'
  if (is.null(X))
  {
    bc <- 0.0
  } else
  {
    bc <- as.vector(X %*% coef(object)$splF)
  }
  sc <- as.vector(B12x %*% coef(object)$splR)
  fit <- mu + bc + sc
  fit
  p.data <- list(x1=x1grid, x2=x2grid)
  L <- list(p.data=p.data, eta=fit, mu=fit)

  return(L)
}


####
####  below, old code, for comparison
####

#' predict for sap2D, without spectral decomposition
#'
#' @inheritParams predict.sap3Dfast
#'
#' @export
predict.sap2Dfast <- function(object,
                              ...,
                              grid) {
  # make predictions on a grid, use min(x1) and max(x1) as in original SAP:
  x1 <- object$x1
  x2 <- object$x2
  x1grid <- seq(min(x1), max(x1), length = grid[1])
  x2grid <- seq(min(x2), max(x2), length = grid[2])

  Bx1 <- spam::as.spam(Bsplines(object$knots1, x1grid))
  Bx2 <- spam::as.spam(Bsplines(object$knots2, x2grid))

  B12x <- Bx1 %x% Bx2

  if (object$scaleX) {
    Xpred <- B12x %*% object$U0
    Xpred[, 1] <- 1.0
  } else {
    X1 <- cbind(1, x1grid)
    X2 <- cbind(1, x2grid)
    Xpred <- X1 %x% X2
  }

  a <- object$a
  p <- ncol(Xpred)
  bc <- as.vector(Xpred %*% a[1:p])
  sc <- as.vector(B12x %*% a[-c(1:p)])
  fit <- bc + sc
  fit
  p.data <- list(x1=x1grid,x2=x2grid)
  L <- list(p.data, eta=fit, mu=fit)
  class(L) <- "predict.sap2Dfast"
  return(L)
}


#' sap2D, without spectral decomposition
#'
#' @inheritParams sap2D
#' @inheritParams LMMsolve
#'
#' @param y ...
#'
#' @export
sap2Dfast <- function(y,
                      x1,
                      x2,
                      nseg,
                      trace = TRUE,
                      tolerance = 1.0e-8,
                      scaleX = FALSE) {
  #nseg <- knots
  n <- length(y)
  s <- proc.time()[3]

  pord <- 2
  degr <- 3
  x1lim <- c(min(x1) - 0.01, max(x1) + 0.01)
  x2lim <- c(min(x2) - 0.01, max(x2) + 0.01)

  knots1 <- PsplinesKnots(x1lim[1], x1lim[2], degree = degr, nseg = nseg[1])
  knots2 <- PsplinesKnots(x2lim[1], x2lim[2], degree = degr, nseg = nseg[2])

  B1 <- spam::as.spam(Bsplines(knots1, x1))
  B2 <- spam::as.spam(Bsplines(knots2, x2))
  q1 <- ncol(B1)
  q2 <- ncol(B2)

  D1 <- spam::diff.spam(diag(q1), diff = pord)
  DtD1 <- crossprod(D1)

  D2 <- spam::diff.spam(diag(q2), diff = pord)
  DtD2 <- crossprod(D2)

  # we have to calculate RowKronecker product only once:
  one.1 <- matrix(1, 1, ncol(B1))
  one.2 <- matrix(1, 1, ncol(B2))
  B12 <- (B1 %x% one.2) * (one.1 %x% B2)

  # calculate the linear/fixed parts:
  U1_null <- cbind(1, scale(1:q1))
  U2_null <- cbind(1, scale(1:q2))

  U1_null <- apply(U1_null, MARGIN = 2, function(x) (x / normVec(x)))
  U2_null <- apply(U2_null, MARGIN = 2, function(x) (x / normVec(x)))

  U_null <- U1_null %x% U2_null
  if (scaleX) {
    X <- B12 %*% U_null
    # take first column as intercept...
    X[,1] <- 1
  } else {
    X1 <- cbind(1, x1)
    X2 <- cbind(1, x2)
    X <- RowKronecker(X1, X2)
  }
  C1 <- spam::spam(x = 0, nrow = q1, ncol = pord)
  C1[1, 1] = C1[q1,2] = 1
  C2 <- spam::spam(x = 0, nrow = q2, ncol = pord)
  C2[1, 1] = C2[q2,2] = 1
  C <- C1 %x% C2

  lRinv <- list()
  lRinv[[1]] = spam::diag.spam(1, n)
  names(lRinv) = "residual"

  CCt <- spam::as.spam(C %*% t(C))
  lGinv <- list()
  lGinv[[1]] <- DtD1 %x% spam::diag.spam(q2) + CCt
  lGinv[[2]] <- spam::diag.spam(q1) %x% DtD2 + CCt

  names(lGinv) <- c("x1", "x2")

  obj <- sparseMixedModels(y = y, X = X, Z = B12, lGinv = lGinv, lRinv = lRinv,
                           maxit = 200, tolerance = tolerance,
                           display = FALSE, trace = trace)

  e <- proc.time()[3]
  cat("Computation time sparse SAP fast:", e - s, "seconds\n")

  L <- list(a = obj$a, edf = obj$ED[-1], knots1 = knots1, knots2 = knots2,
            U0 = U_null, x1 = x1, x2 = x2, scaleX = scaleX)

  class(L) <- "sap2Dfast"
  return(L)
}









