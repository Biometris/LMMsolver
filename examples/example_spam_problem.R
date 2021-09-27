
library(spam)

K <- 100
r <- 2
lambda <- 1.0
B <- diag.spam(r+lambda, nrow=K+1)
B[,1 ] <- r
B[1,] <- r
B[1, 1] <- K*r
B[1:4, 1:4]

cholB1 <- spam::chol(B)
str(cholB1)
cholB2 <- spam::chol(B, memory=list(nnzcolindices=50000))

str(cholB)

sessionInfo()
