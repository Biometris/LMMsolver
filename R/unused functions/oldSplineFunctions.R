####
####  below, old code, for comparison
####

#' predict for spl2D, without spectral decomposition
#'
#' @inheritParams predict.spl3Dfast
#'
#' @export
predict.spl2Dfast <- function(object,
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
  class(L) <- "predict.spl2Dfast"
  return(L)
}


#' spl2D, without spectral decomposition
#'
#' @inheritParams spl2D
#' @inheritParams LMMsolve
#'
#' @param y ...
#'
#' @export
spl2Dfast <- function(y,
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
                           maxit = 200, tolerance = tolerance, trace = trace)

  e <- proc.time()[3]
  cat("Computation time sparse SAP fast:", e - s, "seconds\n")

  L <- list(a = obj$a, edf = obj$ED[-1], knots1 = knots1, knots2 = knots2,
            U0 = U_null, x1 = x1, x2 = x2, scaleX = scaleX)

  class(L) <- "spl2Dfast"
  return(L)
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

  B1 <- Bsplines(knots1, x1)
  B2 <- Bsplines(knots2, x2)
  B3 <- Bsplines(knots3, x3)
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
                           maxit = 200, tolerance = tolerance, trace = trace)

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











