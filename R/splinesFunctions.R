## some help functions, make part of LMM solve library:

# row-wise kronecker product
RowKronecker <-
  function(X1,X2) {
    one.1 <- matrix(1,1,ncol(X1))
    one.2 <- matrix(1,1,ncol(X2))
    kronecker(X1,one.2)*kronecker(one.1,X2)
  }

ndxMatrix <- function(df, lZ, Names)
{
  n <- length(lZ)
  dim <- sapply(lZ,ncol)
  e <- cumsum(dim) + ncol(df)
  s <- e - dim + 1
  lM <- list()
  for (i in 1:n) {
    lM[[i]] <- c(s[i]:e[i])
  }
  names(lM) <- Names
  lM
}

# spectral decomposition of D'D, returns a q x (q-ord) matrix
calcUsc <- function(q, ord)
{
  D <- diff(diag(q), diff=ord)
  DtD <- crossprod(D)
  U <- eigen(DtD)$vectors[,1:(q-ord)]
  d <- eigen(DtD)$values[1:(q-ord)]
  U %*% diag(1/sqrt(d))
}

# equally placed knots:
PsplinesKnots <- function(xmin, xmax, degree, nseg)
{
  dx = (xmax - xmin) / nseg
  knots <- seq(xmin - degree * dx, xmax + degree * dx, by = dx)
  attr(knots, "degree") <- degree
  knots
}

# just to simplify splineDesign
Bsplines <- function(knots, x, deriv = 0)
{
  degree <- attr(knots,"degree")
  B = splineDesign(knots, x, derivs=rep(deriv, length(x)), ord = degree+1)
  B
}
