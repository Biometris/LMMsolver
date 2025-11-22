library(spam)

set.seed(1234)

## tests for RowKronecker

X1 <- spam(x = rnorm(20), nrow = 10, ncol = 2)
X2 <- spam(x = rnorm(40), nrow = 10, ncol = 4)

Z1 <- LMMsolver:::RowKronecker(X1, X2)

X1d <- as.matrix(X1)
X2d <- as.matrix(X2)
Z2 <- LMMsolver:::RowKronecker(X1d, X2d)

expect_equal(as.matrix(Z1), Z2)

X3 <- spam(x = rnorm(24), nrow = 6, ncol = 4)

expect_error(LMMsolver:::RowKronecker(X1 = X1, X2 = X3),
             "X1 and X2 have unequal number of rows.")

## tests for MatrixProduct

A <- spam(x = rnorm(20), nrow = 5, ncol = 4)
B <-spam(x = rnorm(30), nrow = 5, ncol = 6)
expect_error(LMMsolver:::MatrixProduct(A,B),
             "MatrixProduct wrong dimensions")

## tests for GetIntVector:
expect_error(LMMsolver:::GetIntVector(A, "rowpointers", 2),
             "argument ArrayIndex should be 0-based (C/C++) or 1-based (R).",
             fixed=TRUE)

## test for expandGinv
lGinv1 <- diag.spam(1,4)
lGinv2 <- NULL
lGinv3 <- LMMsolver:::expandGinv(lGinv1, lGinv2)
expect_equal(lGinv1, lGinv3)

# test for dlogdet
lC <- list()
lC[[1]] <- diag.spam(1, 3)
ADobj <- LMMsolver:::ADchol(lC)
theta <- c(1 , 2)
expect_error(LMMsolver:::dlogdet(ADobj, theta),
             "wrong length vector theta")

# tests for MatrixProduct
library(Matrix)
A1 <- matrix(data=c(1:12), nrow=3, ncol=4)
A2 <- Matrix(data=c(1:12), nrow=3, ncol=4)
B <- diag.spam(1, 4)
expect_error(LMMsolver:::MatrixProduct(A1, B),
             "Not an S4 object")
expect_error(LMMsolver:::MatrixProduct(A2, B),
             "Both arguments for MatrixProduct should be of class spam")

expect_error(LMMsolver:::chkSplinesFormula("wrongArg"),
             "spline should be a formula of form \"~ spl1D() + ... +  spl1D()\", \"~ spl2D()\" or \"~spl3D()\"",
             fixed=TRUE)


