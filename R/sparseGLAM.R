SparseGLAM <- function(B1x, B2x)
{
  # input B1x, B2x
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
  obj <- list(tG1=tG1_nozeros, tG2=tG2_nozeros, q1=q1, q2=q2, s1=s1, s2=s2, A=A)
  class(obj) <- "SparseGLAM"
  obj
}

calcBtWB <- function(obj, w) {
  z <- LMMsolver:::KronProd2(obj$tG1, obj$tG2, w)
  A <- obj$A
  A@entries <- LMMsolver:::ReArrange(A, obj$q1, obj$q2, obj$s1, obj$s2, z)
  A
}

