#' ---
#' title: Tensor Products in LMMsolver
#' author: Martin Boer, Biometris
#' ---

library(LMMsolver)
suppressPackageStartupMessages(library(spam))
set.seed(1234)

x1min <- 0
x1max <- 7
x2min <- 1
x2max <- 6
x3min <- 7
x3max <- 12
knots1 <- LMMsolver:::PsplinesKnots(x1min, x1max, degree=3, nseg=12)
knots2 <- LMMsolver:::PsplinesKnots(x2min, x2max, degree=3, nseg=17)
knots3 <- LMMsolver:::PsplinesKnots(x3min, x3max, degree=3, nseg=22)

x1 <- seq(x1min, x1max, length=20)
x2 <- seq(x2min, x2max, length=25)
x3 <- seq(x3min, x3max, length=20)
Bx1 <- LMMsolver:::Bsplines(knots1, x1)
Bx2 <- LMMsolver:::Bsplines(knots2, x2)
Bx3 <- LMMsolver:::Bsplines(knots3, x3)

#' Calculation of the form yhat = B %*% a
#' ==================

B <- Bx1 %x% Bx2 %x% Bx3
q <- ncol(Bx1)*ncol(Bx2)*ncol(Bx3)
a <- rnorm(q)

yhat1 <- as.vector(B %*% a)
yhat2 <- LMMsolver:::KronProd(Bx1, Bx2, Bx3, a)
yhat3 <- LMMsolver:::KronProdList(list(Bx1, Bx2, Bx3), a)

all.equal(yhat1, yhat2)
all.equal(yhat1, yhat3)

#' Calculation of the form t(B) %*% y
#' ==================

n <- nrow(Bx1)*nrow(Bx2)*nrow(Bx3)
y <- rnorm(n)

yB1 <- as.vector(t(B) %*% y)
yB2 <- LMMsolver:::KronProd(t(Bx1), t(Bx2), t(Bx3), y)
yB3 <- LMMsolver:::KronProdList(list(t(Bx1), t(Bx2), t(Bx3)), y)

all.equal(yB1, yB2)
all.equal(yB1, yB3)

#' Extended format using RowKronecker
#' ==================

one1 <- rep(1, nrow(Bx1))
one2 <- rep(1, nrow(Bx2))
one3 <- rep(1, nrow(Bx3))
xx1 <- x1 %x% one2 %x% one3
xx2 <- one1 %x% x2 %x% one3
xx3 <- one1 %x% one2 %x% x3

# each Bxx has 4*n elements, total 3*(4*n) = 12n elements
# B has (4^3)*n = 64n elements!
# (not depending on nseg, only on degr=3 of P-splines)
Bxx1 <- LMMsolver:::Bsplines(knots1, xx1)
Bxx2 <- LMMsolver:::Bsplines(knots2, xx2)
Bxx3 <- LMMsolver:::Bsplines(knots3, xx3)

tmp <- LMMsolver:::RowKronecker(Bxx1, Bxx2)
B_ext <- LMMsolver:::RowKronecker(tmp, Bxx3)
all.equal(B, B_ext)


#' Test of ordering data, idea for LMMsolver
#' ==================

df1 <- data.frame(x1=xx1, x2=xx2, x3=xx3)

set.seed(1234)

random_order <- sample(1:n, n)
df2 <- df1[random_order, ]
head(df2)

ord <- order(df2$x1, df2$x2, df2$x3)
df3 <- df2[ord, ]

all.equal(df1, df3)

# This could be check for grid in LMMsolver,
# after ordering.
all.equal(df3$x1, as.vector(xx1))
all.equal(df3$x2, as.vector(xx2))
all.equal(df3$x3, as.vector(xx3))
