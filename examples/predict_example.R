#' ---
#' title: Predict using sparse inverse
#' author: Martin Boer, Biometris
#' ---

library(LMMsolver)
library(JOPS)
suppressPackageStartupMessages(library(spam))
set.seed(1234)

#' One dimension, see section 2.1 Currie 2006
#' ==================

xmin <- 0
xmax <- 7
knots <- LMMsolver:::PsplinesKnots(xmin, xmax, degree=3, nseg=5)

n <- 100
x <- seq(xmin, xmax, length=n)
B <- LMMsolver:::Bsplines(knots, x)
q <- ncol(B)
w <- runif(n)
W <- diag(w)

BtWB <- t(B) %*% W %*% B
dim(BtWB)

cholC <- chol(BtWB)
LMMsolver:::PrintCholesky(cholC)
display(t(cholC))

Cinv <- solve(BtWB)

D1 <- as.spam(matrix(data=c(1,2,rep(0,6)),ncol=8,nrow=1))
D2 <- as.spam(matrix(data=c(0,1,rep(0,4),2,3),ncol=8,nrow=1))
D3 <- as.spam(matrix(data=c(2,0,2,0,0,0,0,3),ncol=8,nrow=1))
D4 <- as.spam(matrix(data=c(0,1,0,0,0,0,0,1),ncol=8,nrow=1))
D5 <- as.spam(matrix(data=c(1,1,0,0,0,0,0,1),ncol=8,nrow=1))
D <- rbind(D1, D2, D3, D4, D5)

updateH <- function(tDp, i, j)
{
  d <- ncol(tDp)
  z <- rep(0,d)
  s <- tDp@rowpointers[i]
  e <- tDp@rowpointers[i+1] - 1
  if (s>e) return(z)
  rowndx1 <- tDp@colindices[c(s:e)]
  val1 <- tDp@entries[c(s:e)]
  s <- tDp@rowpointers[j]
  e <- tDp@rowpointers[j+1] - 1
  if (s>e) return(z)
  rowndx2 <- tDp@colindices[c(s:e)]
  val2 <- tDp@entries[c(s:e)]

  col <- intersect(rowndx1, rowndx2)
  u <- val1[which(rowndx1 %in% col)]*val2[which(rowndx2 %in% col)]
  d <- ncol(tDp)
  z[col] <- u
  z
}

# make predictions on a grid:
x0 <- seq(xmin,xmax,length=20)
Bx0 <- LMMsolver:::Bsplines(knots, x0)

# important to add Bx0!!
C = BtWB + 0 * crossprod(Bx0)
cholC <- chol(C)
cholC@entries <- LMMsolver:::partialDerivCholesky(cholC)
A <- spam::as.spam(cholC)

# rorder B
p <- cholC@pivot
Bp <- Bx0[, p]
display(Bp)
d <- nrow(Bx0)
# by foot
tBp <- t(Bp)
s <- rep(0, d)
q <- ncol(BtWB)
for (i in 1:q) {
  for (j in 1:q) {
    alpha <- A[i,j]
    z <- updateH(tBp, i, j)
    s <- s + alpha*z
  }
}
s

s2 <- diag(Bx0 %*% Cinv %*% t(Bx0))
s3 <- LMMsolver:::diagXCinvXt(chol(C), tBp)
range(s-s2)
range(s-s3)

se1 <- LMMsolver:::calcStandardErrors(C, Bx0, NewMethod=FALSE)
se2 <- LMMsolver:::calcStandardErrors(C, Bx0, NewMethod=TRUE)
range(s-se1^2)
range(s-se2^2)
