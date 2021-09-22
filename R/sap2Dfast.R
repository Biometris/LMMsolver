#' sap2D, without spectral decomposition
#'
#' sap2D, without spectral decomposition
#'
#' @param x1 ...
#' @param x2 ...
#' @param knots ...
#' @param scaleX ...
#'
#' @export
sap2D <- function(x1,
                  x2,
                  knots,
                  scaleX = FALSE) {
  nseg <- knots

  pord <- 2
  degr <- 3
  #x1lim <- c(min(x1) - 0.01, max(x1) + 0.01)
  #x2lim <- c(min(x2) - 0.01, max(x2) + 0.01)
  x1lim <- c(min(x1), max(x1))
  x2lim <- c(min(x2), max(x2))

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

  if (scaleX) {
    ## calculate the linear/fixed parts.
    U1_null <- cbind(1, scale(1:q1))
    U2_null <- cbind(1, scale(1:q2))

    U1_null <- apply(U1_null, MARGIN = 2, function(x) (x / normVec(x)))
    U2_null <- apply(U2_null, MARGIN = 2, function(x) (x / normVec(x)))

    U_null <- U1_null %x% U2_null

    X <- B12 %*% U_null
    ## Remove intercept column to avoid singularity problems.
    X <- X[, -1]
  } else {
    X1 <- cbind(1, x1)
    X2 <- cbind(1, x2)
    X <- RowKronecker(X1, X2)
    ## Remove intercept column to avoid singularity problems.
    X <- X[, -1]
  }
  C1 <- spam::spam(x = 0, nrow = q1, ncol = pord)
  C1[1, 1] = C1[q1,2] = 1
  C2 <- spam::spam(x = 0, nrow = q2, ncol = pord)
  C2[1, 1] = C2[q2,2] = 1
  C <- C1 %x% C2

  CCt <- spam::tcrossprod(C)

  lGinv <- list()
  lGinv[[1]] <- DtD1 %x% spam::diag.spam(q2) + CCt
  lGinv[[2]] <- spam::diag.spam(q1) %x% DtD2 + CCt

  names(lGinv) <- c("x1", "x2")

  return(list(X = X, Z = B12, lGinv = lGinv))
}

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
                      knots,
                      trace = TRUE,
                      tolerance = 1.0e-8,
                      scaleX = FALSE) {
  nseg <- knots
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
