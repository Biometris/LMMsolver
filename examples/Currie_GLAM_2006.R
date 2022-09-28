#' ---
#' title: Some analysis Currie 2006 paper
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

G <- LMMsolver:::RowKronecker(B, B)

# see equation (2.3) and the one below, p. 262 Currie et al:
z1 <- as.vector(as.matrix(BtWB))
z2 <- as.vector(t(G) %*% w)
all.equal(z1, z2)

# example by foot, to understand extension to higher dimensions:
i <- 3
j <- 4
ndx <- (j-1)*q + (i-1) + 1
z2[ndx]
BtWB[i, j]
as.numeric(t(B[, i]) %*% W %*% B[, j])

# tG is sparse, with empty rows...
tG <- t(G)
s <- which(diff(tG@rowpointers)>0)
tG2 <- tG[s, ]
tG2[s==ndx] %*% w

# w=1, just to have the structure of BtWB
BtWB2 <- crossprod(B)
BtWB2@entries <- as.vector(tG2 %*% w)
all.equal(BtWB, BtWB2)

Cinv <- solve(BtWB)

# calculate the diagonal elements B (BtWB)^{-1} B'
cholC <- chol(BtWB)
## calculate the partial derivatives of Cholesky
cholC@entries <- LMMsolver:::partialDerivCholesky(cholC)
A <- spam::as.spam(cholC)
A <- A[cholC@invpivot, cholC@invpivot]
a <- as.vector(as.matrix(A))

dg1 <- diag(B %*% Cinv %*% t(B))
dg2 <- as.vector(G %*% as.vector(Cinv))
dg3 <- as.vector(G %*% a)
all.equal(dg1, dg2)
all.equal(dg1, dg3)

#' Two dimensions, see section 2.1 Currie 2006
#' ==================

x1min <- 0
x1max <- 7
x2min <- 3
x2max <- 12
knots1 <- LMMsolver:::PsplinesKnots(x1min, x1max, degree=3, nseg=7)
knots2 <- LMMsolver:::PsplinesKnots(x2min, x2max, degree=3, nseg=12)

n1 <- 50
n2 <- 60
n <- n1*n2
x1 <- seq(x1min, x1max, length=n1)
x2 <- seq(x2min, x2max, length=n2)

B1x <- LMMsolver:::Bsplines(knots1, x1)
B2x <- LMMsolver:::Bsplines(knots2, x2)
q1 <- ncol(B1x)
q2 <- ncol(B2x)
B <- B1x %x% B2x
q <- ncol(B)

w <- runif(n)
W <- diag.spam(w)

BtWB <- t(B) %*% W %*% B
dim(BtWB)

# use array trick:
G1 <- LMMsolver:::RowKronecker(B1x, B1x)
G2 <- LMMsolver:::RowKronecker(B2x, B2x)
tG1 <- t(G1)
tG2 <- t(G2)
Wa <- matrix(data=w, n2, n1)
Q <- as.matrix(tG2 %*% Wa %*% G1)
Q_org <- Q
dim(Q) <- c(q2, q2, q1, q1)
Q <- aperm(Q, c(1, 3, 2, 4))
dim(Q) <- c(q1 * q2, q1 * q2)
all.equal(as.matrix(BtWB), Q)

# use C++ version, without reordering
z3 <- LMMsolver:::KronProd2(tG1, tG2, w)
Q_org2 <- matrix(data=z3, q2^2, q1^2)
all.equal(Q_org2, Q_org)

#' Example by foot
#' ==================

i1 <- 2
j1 <- 3
i2 <- 3
j2 <- 4

tB1x <- t(B1x)
tB2x <- t(B2x)

# method 1 to get value of BtWB[i,j]
i <- (i1-1)*q2 + i2
j <- (j1-1)*q2 + j2
BtWB[i, j]

# method 2 to get value of BtWB[i,j]
A1 <- as.matrix((tB1x[i1, ] %x% tB2x[i2, ]))
A2 <- as.matrix((tB1x[j1, ] %x% tB2x[j2, ]))
as.numeric(A1 %*% W %*% t(A2))

# method 3: more efficient method to get BtWB[i,j]
ndx1 <- (j1-1)*q1 + i1
ndx2 <- (j2-1)*q2 + i2
tG <- tG1[ndx1, ] %x% tG2[ndx2, ]
as.numeric(tG %*% w)

# method 4: use C++/Rpp
z <- LMMsolver:::KronProd2(tG1, tG2, w)
#ndx <- (i1-1)*(q1*q2^2) + (j1-1)*(q2^2) + (i2-1)*q2 + (j2-1) + 1
ndx <- ((i1-1)*q1 + (j1-1))*(q2^2) + ((i2-1)*q2 + (j2-1)) + 1
z[ndx]

#' Sparse solution
#' ==================

# make the SparseGLAM object:
W <- diag(w)

# Prpare the penalty matrices
Dx1 = diff(diag(q1), diff = 2)
Dx2 = diff(diag(q2), diff = 2)
lambdax = lambday = 1
Px1 = lambdax * t(Dx1) %*% Dx1
Px2 = lambday * t(Dx2) %*% Dx2
P = kronecker(Px2, diag(q1)) + kronecker(diag(q2), Px1)
P <- as.spam(P)

C <- t(B) %*% W %*% B + P
obj <- LMMsolver:::sparseGLAM(B1x, B2x)
# calculate BtWB using weight w:
M <- LMMsolver:::calcBtWB(obj, w)
all.equal(BtWB, M)

Cinv <- solve(C)
## Equivalent to diag(B %*% Cinv %*% t(B))
dg1 <- spam::rowSums.spam((B %*% Cinv) * B)
dg2 <- LMMsolver:::calcB_Cinv_Bt(obj, C)
all.equal(dg1, dg2)
