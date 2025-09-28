library(spam)

set.seed(1234)

X1 <- spam(x = rnorm(20), nrow = 10, ncol = 2)
X2 <- spam(x = rnorm(40), nrow = 10, ncol = 4)

Z1 <- LMMsolver:::RowKronecker(X1, X2)

X1d <- as.matrix(X1)
X2d <- as.matrix(X2)
Z2 <- LMMsolver:::RowKronecker(X1d, X2d)

expect_equal(as.matrix(Z1), Z2)

