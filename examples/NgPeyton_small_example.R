library(spam)
library(LMMsolver)

# nice small example with four supernodes:
n1 = 2
n2 = 5
D1 = diff(diag(n1),diff=1)
D2 = diff(diag(n2),diff=1)

A1 = diag(n1) + t(D1)%*%D1
A2 = diag(n2) + t(D2)%*%D2
B = kronecker(A1, A2)
B
B = as.spam(B)
U = chol(B)
U@supernodes

diff(U@supernodes)

# use lower triangle:
L <- t(U)
colpointers <- U@rowpointers
rowpointers <- U@colpointers
rowindices <- U@colindices

U@supernodes
colpointers
rowpointers
rowindices

display(L)
nNodes = length(U@supernodes)-2
for (i in 1:nNodes)
  abline(v=U@supernodes[i+1]-0.5, col='red')

lP <- list()
lP[[1]] <- B

obj0 <- LMMsolver:::ADchol(lP)
obj1 <- LMMsolver:::ADcholNgPeyton(lP)

lambda = 1.0
det0 <- LMMsolver:::logdet(obj0, lambda=lambda)
det1 <- LMMsolver:::logdetNgPeyton(obj1, lambda=lambda)
det2 <- as.numeric(determinant(lambda*B)$modulus)

det0
det1
det2
det0-det2
det1-det2

# output outline new version:
LMMsolver:::dlogdetNgPeyton(obj1,lambda=lambda)

