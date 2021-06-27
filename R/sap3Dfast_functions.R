#' sap3D, without spectral decomposition
#'
#' @export
sap3Dfast <- function(y, x1, x2, x3, knots, trace=TRUE, thr=1.0e-8, scaleX=FALSE)
{
  nseg <- knots
  s <- proc.time()[3]

  pord <- 2
  degr <- 3
  x1lim <- c(min(x1)-0.01,max(x1)+0.01)
  x2lim <- c(min(x2)-0.01,max(x2)+0.01)
  x3lim <- c(min(x3)-0.01,max(x3)+0.01)

  knots1 <- PsplinesKnots(x1lim[1], x1lim[2], degree=degr, nseg=nseg[1])
  knots2 <- PsplinesKnots(x2lim[1], x2lim[2], degree=degr, nseg=nseg[2])
  knots3 <- PsplinesKnots(x3lim[1], x3lim[2], degree=degr, nseg=nseg[3])

  B1 <- as.spam(Bsplines(knots1, x1))
  B2 <- as.spam(Bsplines(knots2, x2))
  B3 <- as.spam(Bsplines(knots3, x3))
  q1 <- ncol(B1)
  q2 <- ncol(B2)
  q3 <- ncol(B3)

  D1 <- diff.spam(diag(q1), diff=pord)
  DtD1 <- crossprod(D1)

  D2 <- diff.spam(diag(q2), diff=pord)
  DtD2 <- crossprod(D2)

  D3 <- diff.spam(diag(q3), diff=pord)
  DtD3 <- crossprod(D3)

  # we have to calculate RowKronecker product only once:
  one.1 <- matrix(1, 1, ncol(B1))
  one.2 <- matrix(1, 1, ncol(B2))
  one.3 <- matrix(1, 1, ncol(B3))
  B123 <- (B1 %x% one.2 %x% one.3) * (one.1 %x% B2 %x% one.3) * (one.1 %x% one.2 %x% B3)
  summary(B123)

  # calculate the linear/fixed parts:
  U1_null <- cbind(1, scale(1:q1))
  U2_null <- cbind(1, scale(1:q2))
  U3_null <- cbind(1, scale(1:q3))

  norm_vec <- function(x) { return(sqrt(sum(x^2)))}
  U1_null <- apply(U1_null, MARGIN=2, function(x) (x/norm_vec(x)))
  U2_null <- apply(U2_null, MARGIN=2, function(x) (x/norm_vec(x)))
  U3_null <- apply(U3_null, MARGIN=2, function(x) (x/norm_vec(x)))

  U_null <- U1_null %x% U2_null %x% U3_null
  if (scaleX) {
    X <- B123 %*% U_null
    # take first column as intercept...
    X[,1] <- 1
  } else {
    X1 <- cbind(1, x1)
    X2 <- cbind(1, x2)
    X3 <- cbind(1, x3)
    X = RowKronecker(RowKronecker(X1,X2), X3)
  }
  C1 <- spam(x=0, nrow=q1, ncol=pord)
  C1[1,1] = C1[q1,2] = 1
  C2 <- spam(x=0, nrow=q2, ncol=pord)
  C2[1,1] = C2[q2,2] = 1
  C3 <- spam(x=0, nrow=q3, ncol=pord)
  C3[1,1] = C3[q3,2] = 1
  C <- C1 %x% C2 %x% C3

  lRinv <- list()
  lRinv[[1]] = diag.spam(1, n)
  names(lRinv) = "residual"

  CCt <- as.spam(C %*% t(C))
  lGinv <- list()
  lGinv[[1]] <- DtD1 %x% diag.spam(q2) %x% diag.spam(q3) + CCt
  lGinv[[2]] <- diag.spam(q1) %x% DtD2 %x% diag.spam(q3) + CCt
  lGinv[[3]] <- diag.spam(q1) %x% diag.spam(q2) %x% DtD3 + CCt

  names(lGinv) <- c('x1','x2','x3')

  obj = sparseMixedModels(y, X, B123, lGinv, lRinv,
                          maxiter=200, eps=thr, display=FALSE, monitor=trace)

  e <- proc.time()[3]
  cat("Computation time sparse SAP fast:", e-s, "seconds\n")

  L = list(a=obj$a, edf=obj$ED[-1], knots1=knots1, knots2=knots2,knots3=knots3,
           U0=U_null, x1=x1,x2=x2,x3=x3,scaleX=scaleX)

  class(L) = "sap3Dfast"
  L
}

#' predict for sap3D, without spectral decomposition
#'
#' @export
predict.sap3Dfast <- function(object, grid)
{
  # make predictions on a grid, use min(x1) and max(x1) as in original SAP:
  x1 <- object$x1
  x2 <- object$x2
  x3 <- object$x3
  x1grid <- seq(min(x1), max(x1), length=grid[1])
  x2grid <- seq(min(x2), max(x2), length=grid[2])
  x3grid <- seq(min(x3), max(x3), length=grid[3])

  Bx1 <- as.spam(Bsplines(object$knots1, x1grid))
  Bx2 <- as.spam(Bsplines(object$knots2, x2grid))
  Bx3 <- as.spam(Bsplines(object$knots3, x3grid))

  B123x <- Bx1 %x% Bx2 %x% Bx3

  if (object$scaleX) {
    Xpred <- B123x %*% object$U0
    Xpred[,1] <- 1.0
  } else {
    X1 <- cbind(1, x1grid)
    X2 <- cbind(1, x2grid)
    X3 <- cbind(1, x3grid)
    Xpred <- X1 %x% X2 %x% X3
  }

  a <- object$a
  p <- ncol(Xpred)
  bc <- as.vector(Xpred %*% a[1:p])
  sc <- as.vector(B123x %*% a[-c(1:p)])
  fit <- bc + sc
  fit
  p.data <- list(x1=x1grid,x2=x2grid,x3=x3grid)
  L = list(p.data, eta=fit, mu=fit)
  class(L) = "predict.sap3Dfast"
  L
}

