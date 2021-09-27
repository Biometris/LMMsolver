rm(list=ls())
library(spam)

K <- 10
r <- 2
lambda <- 1.0
B <- diag.spam(r+lambda, nrow=K+1)
B[,1 ] <- r
B[1,] <- r
B[1, 1] <- K*r
B[1:4, 1:4]

# for K <= 500 ok, for K > 500 the function crashes often,
# there seems to be problem with memory allocation:
cholB1 <- spam::chol(B)
str(cholB1)

# this is ok...
cholB2 <- spam::chol(B, memory=list(nnzcolindices=50000))

str(cholB2)

sessionInfo()
