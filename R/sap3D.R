#' fit 3D P-splines
#'
#' fit 3D P-splines using sparse implementation.
#'
#' @param x1 numerical vector containing the values of \code{x1} covariate.
#' @param x2 numerical vector containing the values of \code{x2} covariate.
#' @param x3 numerical vector containing the values of \code{x3} covariate.

#' @param nseg number of segments
#' @param scaleX logical, scale fixed effect or not. Default is TRUE,
#' no scaling.
#' @param pord order of penalty, default \code{pord=2}
#' @param degree degree of B-spline basis, default \code{degree=3}
#' @param x1lim numerical vector of length 2 containing the domain of covariate
#' \code{x1} where the knots should be placed. Default set to \code{NULL} (covariate range).
#' @param x2lim numerical vector of length 2 containing the domain of covariate
#' \code{x2} where the knots should be placed. Default set to \code{NULL} (covariate range).
#' @param x3lim numerical vector of length 2 containing the domain of covariate
#' \code{x3} where the knots should be placed. Default set to \code{NULL} (covariate range).
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
spl3D <- function(x1,
                  x2,
                  x3,
                  nseg,
                  pord = 2,
                  degree = 3,
                  scaleX = TRUE,
                  x1lim = NULL,
                  x2lim = NULL,
                  x3lim = NULL) {
  if (is.null(x1lim))
    x1lim <- c(min(x1), max(x1))
  if (is.null(x2lim))
    x2lim <- c(min(x2), max(x2))
  if (is.null(x3lim))
    x3lim <- c(min(x3), max(x3))

  knots <- list(1)
  knots[[1]] <- PsplinesKnots(x1lim[1], x1lim[2], degree = degree, nseg = nseg[1])
  knots[[2]] <- PsplinesKnots(x2lim[1], x2lim[2], degree = degree, nseg = nseg[2])
  knots[[3]] <- PsplinesKnots(x3lim[1], x3lim[2], degree = degree, nseg = nseg[3])

  B1 <- spam::as.spam(Bsplines(knots[[1]], x1))
  B2 <- spam::as.spam(Bsplines(knots[[2]], x2))
  B3 <- spam::as.spam(Bsplines(knots[[3]], x3))

  q1 <- ncol(B1)
  q2 <- ncol(B2)
  q3 <- ncol(B3)

  # we have to calculate RowKronecker product only once:
  one.1 <- matrix(1, 1, ncol(B1))
  one.2 <- matrix(1, 1, ncol(B2))
  one.3 <- matrix(1, 1, ncol(B3))
  B123 <- (B1 %x% one.2 %x% one.3) *
          (one.1 %x% B2 %x% one.3) *
          (one.1 %x% one.2 %x% B3)

  DtD1 <- constructPenalty(q1, pord)
  DtD2 <- constructPenalty(q2, pord)
  DtD3 <- constructPenalty(q3, pord)

  X1 <- constructX(B1, x1, scaleX, pord)
  X2 <- constructX(B2, x2, scaleX, pord)
  X3 <- constructX(B3, x3, scaleX, pord)

  X <- RowKronecker(RowKronecker(X1, X2), X3)

  ## Remove intercept column to avoid singularity problems.
  X  <- removeIntercept(X)

  CCt1 <- constructCCt(q1, pord)
  CCt2 <- constructCCt(q2, pord)
  CCt3 <- constructCCt(q3, pord)

  CCt <- CCt1 %x% CCt2 %x% CCt3

  lGinv <- list()
  lGinv[[1]] <- DtD1 %x% spam::diag.spam(q2) %x% spam::diag.spam(q3) + CCt
  lGinv[[2]] <- spam::diag.spam(q1) %x% DtD2 %x% spam::diag.spam(q3) + CCt
  lGinv[[3]] <- spam::diag.spam(q1) %x% spam::diag.spam(q2) %x% DtD3 + CCt

  names(lGinv) <- c("s(x1)", "s(x2)", "s(x3)")

  if (is.null(X))
  {
    dim.f = NULL
    term.labels.f = NULL
  } else {
    dim.f = c(ncol(X))
    term.labels.f = c('splF')
  }
  dim.r = c(ncol(B123))
  term.labels.r = c('splR')
  return(list(X = X, Z = B123, lGinv = lGinv, knots = knots,
              dim.f=dim.f, dim.r=dim.r, term.labels.f=term.labels.f,
              term.labels.r=term.labels.r,x1=x1,x2=x2,x3=x3,
              pord=pord, degree=degree, scaleX=scaleX))
}



##
## old code below (check for making predictions)
##

#' spl3D, without spectral decomposition
#'
#' @inheritParams spl3D
#' @inheritParams LMMsolve
#'
#' @param y ...
#'
#' @export
spl3Dfast <- function(y,
                      x1,
                      x2,
                      x3,
                      nseg,
                      trace = TRUE,
                      tolerance = 1.0e-8,
                      scaleX = TRUE) {
  #nseg <- knots
  n <- length(y)
  s <- proc.time()[3]

  pord <- 2
  degr <- 3
  x1lim <- c(min(x1) - 0.01, max(x1) + 0.01)
  x2lim <- c(min(x2) - 0.01, max(x2) + 0.01)
  x3lim <- c(min(x3) - 0.01, max(x3) + 0.01)

  knots1 <- PsplinesKnots(x1lim[1], x1lim[2], degree = degr, nseg = nseg[1])
  knots2 <- PsplinesKnots(x2lim[1], x2lim[2], degree = degr, nseg = nseg[2])
  knots3 <- PsplinesKnots(x3lim[1], x3lim[2], degree = degr, nseg = nseg[3])

  B1 <- spam::as.spam(Bsplines(knots1, x1))
  B2 <- spam::as.spam(Bsplines(knots2, x2))
  B3 <- spam::as.spam(Bsplines(knots3, x3))
  q1 <- ncol(B1)
  q2 <- ncol(B2)
  q3 <- ncol(B3)

  D1 <- spam::diff.spam(diag(q1), diff = pord)
  DtD1 <- crossprod(D1)

  D2 <- spam::diff.spam(diag(q2), diff = pord)
  DtD2 <- crossprod(D2)

  D3 <- spam::diff.spam(diag(q3), diff = pord)
  DtD3 <- crossprod(D3)

  # we have to calculate RowKronecker product only once:
  one.1 <- matrix(1, 1, ncol(B1))
  one.2 <- matrix(1, 1, ncol(B2))
  one.3 <- matrix(1, 1, ncol(B3))
  B123 <- (B1 %x% one.2 %x% one.3) *
    (one.1 %x% B2 %x% one.3) *
    (one.1 %x% one.2 %x% B3)
  summary(B123)

  # calculate the linear/fixed parts:
  U1_null <- cbind(1, scale(1:q1))
  U2_null <- cbind(1, scale(1:q2))
  U3_null <- cbind(1, scale(1:q3))

  U1_null <- apply(U1_null, MARGIN = 2, function(x) (x / normVec(x)))
  U2_null <- apply(U2_null, MARGIN = 2, function(x) (x / normVec(x)))
  U3_null <- apply(U3_null, MARGIN = 2, function(x) (x / normVec(x)))

  U_null <- U1_null %x% U2_null %x% U3_null
  if (scaleX) {
    X <- B123 %*% U_null
    # take first column as intercept...
    X[, 1] <- 1
  } else {
    X1 <- cbind(1, x1)
    X2 <- cbind(1, x2)
    X3 <- cbind(1, x3)
    X = RowKronecker(RowKronecker(X1, X2), X3)
  }
  C1 <- spam::spam(x = 0, nrow = q1, ncol = pord)
  C1[1, 1] <- C1[q1, 2] <- 1
  C2 <- spam::spam(x = 0, nrow = q2, ncol = pord)
  C2[1, 1] <- C2[q2, 2] <- 1
  C3 <- spam::spam(x = 0, nrow = q3, ncol = pord)
  C3[1, 1] <- C3[q3, 2] <- 1
  C <- C1 %x% C2 %x% C3

  lRinv <- list()
  lRinv[[1]] = spam::diag.spam(1, n)
  names(lRinv) = "residual"

  CCt <- spam::as.spam(C %*% t(C))
  lGinv <- list()
  lGinv[[1]] <- DtD1 %x% spam::diag.spam(q2) %x% spam::diag.spam(q3) + CCt
  lGinv[[2]] <- spam::diag.spam(q1) %x% DtD2 %x% spam::diag.spam(q3) + CCt
  lGinv[[3]] <- spam::diag.spam(q1) %x% spam::diag.spam(q2) %x% DtD3 + CCt

  names(lGinv) <- c("x1", "x2", "x3")

  obj <- sparseMixedModels(y, X, B123, lGinv, lRinv,
                           maxit = 200, tolerance = tolerance, display = FALSE,
                           trace = trace)

  e <- proc.time()[3]
  cat("Computation time sparse SAP fast:", e - s, "seconds\n")

  L <- list(a = obj$a, edf = obj$ED[-1], knots1 = knots1, knots2 = knots2,
            knots3 = knots3, U0 = U_null, x1 = x1, x2 = x2, x3 = x3,
            scaleX = scaleX)

  class(L) <- "spl3Dfast"
  return(L)
}

#' predict for spl3D, without spectral decomposition
#'
#' @param object ...
#' @param ... ...
#' @param newdata ...
#' @param grid ...
#'
#' @export
predict.spl3Dfast <- function(object,
                              ...,
                              newdata,
                              grid = c(30, 30, 30)) {
  if (!missing(newdata)) {
    if (!is.data.frame(newdata)) {
      stop("newdata should be a dataframe")
    }
    if (!all(c("x1", "x2", "x3") %in% names(newdata))) {
      stop("Not all needed variables (x1,x2,x3) are supplied in newdata")
    }

    x1p <- newdata$x1
    x2p <- newdata$x2
    x3p <- newdata$x3

    Bx1 <- spam::as.spam(Bsplines(object$knots1, x1p))
    Bx2 <- spam::as.spam(Bsplines(object$knots2, x2p))
    Bx3 <- spam::as.spam(Bsplines(object$knots3, x3p))
    B123x <- RowKronecker(RowKronecker(Bx1, Bx2), Bx3)

  } else { # grid....
    # make predictions on a grid, use min(x1) and max(x1) as in original SAP:
    x1p <- seq(min(object$x1), max(object$x1), length = grid[1])
    x2p <- seq(min(object$x2), max(object$x2), length = grid[2])
    x3p <- seq(min(object$x3), max(object$x3), length = grid[3])
    Bx1 <- spam::as.spam(Bsplines(object$knots1, x1p))
    Bx2 <- spam::as.spam(Bsplines(object$knots2, x2p))
    Bx3 <- spam::as.spam(Bsplines(object$knots3, x3p))
    B123x <- Bx1 %x% Bx2 %x% Bx3
  }


  if (object$scaleX) {
    Xpred <- B123x %*% object$U0
    Xpred[, 1] <- 1.0
  } else {
    X1 <- cbind(1, x1p)
    X2 <- cbind(1, x2p)
    X3 <- cbind(1, x3p)
    if (!missing(newdata)) {
      Xpred <- RowKronecker(RowKronecker(X1, X2), X3)
    } else {
      Xpred <- X1 %x% X2 %x% X3
    }
  }

  a <- object$a
  p <- ncol(Xpred)
  bc <- as.vector(Xpred %*% a[1:p])
  sc <- as.vector(B123x %*% a[-c(1:p)])
  fit <- bc + sc

  if (missing(newdata)) {
    p.data <- list(x1 = x1p, x2 = x2p, x3 = x3p)
  } else {
    p.data <- newdata
  }
  L <- list(p.data = p.data, eta = fit, mu = fit)
  class(L) <- "predict.spl3Dfast"
  return(L)
}
