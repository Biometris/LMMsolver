
sparseGLAM <- function(B1x, B2x)
{
  q1 <- ncol(B1x)
  q2 <- ncol(B2x)
  G1 <- LMMsolver:::RowKronecker(B1x, B1x)
  G2 <- LMMsolver:::RowKronecker(B2x, B2x)
  tG1 <- t(G1)
  tG2 <- t(G2)
  s1 <- which(diff(tG1@rowpointers) != 0)
  s2 <- which(diff(tG2@rowpointers) != 0)
  tG1_nozeros <- tG1[s1, ]
  tG2_nozeros <- tG2[s2, ]
  A <- spam::crossprod.spam(B1x) %x% spam::crossprod.spam(B2x)
  ord <- LMMsolver:::getOrder(A, q1, q2, s1, s2)

  obj <- list(B1x=B1x, B2x=B2x,G1=G1, G2=G2,
              tG1=tG1_nozeros, tG2=tG2_nozeros,
              q1=q1, q2=q2, s1=s1, s2=s2, A=A,ord=ord)
  class(obj) <- "sparseGLAM"
  obj
}

calcBtWB <- function(obj, w) {
  z <- LMMsolver:::KronProd2(obj$tG1, obj$tG2, w)
  A <- obj$A
  A@entries <- z[obj$ord]
  A
}

calcBtY <- function(obj, y) {
  z <- LMMsolver:::KronProd2(t(obj$B1x), t(obj$B2x), y)
  z
}

calcBa <- function(obj, a) {
  z <- LMMsolver:::KronProd2(obj$B1x, obj$B2x, a)
  z
}

calcB_Cinv_Bt <- function(obj, C, w=1) {
  cholC <- chol(C)
  q1 <- obj$q1
  q2 <- obj$q2

  ## calculate the partial derivatives of Cholesky
  cholC@entries <- LMMsolver:::partialDerivCholesky(cholC)

  A <- spam::as.spam(cholC)
  A <- A[cholC@invpivot, cholC@invpivot]
  A <- as.matrix(A)
  dim(A) <- c(q2, q1, q2, q1)
  A <- aperm(A, c(1, 3, 2, 4))
  dim(A) <- c(q2 * q2, q1 * q1)
  a <- as.vector(A)

  x <- LMMsolver:::KronProd2(obj$G1, obj$G2, a)*w
  x
}



