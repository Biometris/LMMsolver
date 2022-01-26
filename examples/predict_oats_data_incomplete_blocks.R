# predict example...
rm(list=ls())

library(agridat)
library(LMMsolver)
library(ggplot2)
library(asreml)
library(spam)

# calculate trace of matrix M:
tr <- function(M) { sum(diag(M)) }

data(john.alpha)
dat <- john.alpha
# remove one observation
# dat <- dat[-20,]

obj1 <- asreml(fixed = yield~rep+gen,
                 random=~rep:block,
                 data = dat)

obj2 <- LMMsolve(fixed = yield~rep+gen,
                random=~rep:block,
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
summary(obj2)

# asreml
coef(obj1, list=TRUE)$'rep:block'
# labels missing
coef(obj2)$'rep:block'

displayMME(obj2)

obj1$loglik
obj2$logL

# predictions using asreml
pred1 <- predict(obj1, classify='gen',sed=TRUE)
head(pred1$pvals[1:4,])
pred1$sed[1:3,1:3]

#
# make predictions by foot ...
#
N <- nrow(dat)
y <- dat$yield
X <- model.matrix(~rep+gen, dat)
Z <- model.matrix(~-1+rep:block,dat)
p <- ncol(X)
q <- ncol(Z)
W <- cbind(X,Z)
Rinv <- diag.spam(1, N)
Ginv <- diag.spam(1, q)
P <- list()
P[[1]] <- diag.spam(c(rep(0, p),rep(1, q)))
P[[2]] <- as.spam(t(W) %*% W)
sapply(P,FUN=dim)

# recalculate u, should be saved in LMMsolve object:
C <- obj2$C
theta <- obj2$theta
WtY <- theta[2] * t(W) %*% y
u <- solve(C, WtY)

# Replace by sparse inverse...
Cinv <- solve(C)

Nrep <- nlevels(dat$rep)
Ngen <- nlevels(dat$gen)
Nblk <- nlevels(dat$block)*Nrep

# predictions plus standard errors,
# averaging over fixed rep, ignoring random rep:block
# to make this general, we have to store fixed/random info in LMMsolve object
Dg1 <- matrix(c(1, rep(1/Nrep, Nrep-1), 0, 0, rep(0, Ngen-3), rep(0, Nblk)), nrow=1)
Dg2 <- matrix(c(1, rep(1/Nrep, Nrep-1), 1, 0, rep(0, Ngen-3), rep(0, Nblk)), nrow=1)
Dg3 <- matrix(c(1, rep(1/Nrep, Nrep-1), 0, 1, rep(0, Ngen-3), rep(0, Nblk)), nrow=1)

p1 <- Dg1 %*% u
p2 <- Dg2 %*% u
p3 <- Dg3 %*% u
se1 <- sqrt(tr(Cinv %*% t(Dg1) %*% Dg1))
se2 <- sqrt(tr(Cinv %*% t(Dg2) %*% Dg2))
se3 <- sqrt(tr(Cinv %*% t(Dg3) %*% Dg3))

pred2 <- c(p1, p2, p3)
pred2se <- c(se1, se2, se3)

# sed
Dg12 <- Dg1-Dg2
Dg13 <- Dg1-Dg3
Dg23 <- Dg2-Dg3
A <- matrix(data=NA,ncol=3,nrow=3)
A[1,2] <- A[2,1] <- sqrt(tr(Cinv %*% t(Dg12) %*% Dg12))
A[1,3] <- A[3,1] <- sqrt(tr(Cinv %*% t(Dg13) %*% Dg13))
A[2,3] <- A[3,2] <- sqrt(tr(Cinv %*% t(Dg23) %*% Dg23))
A

# compare results with asreml
pred2 - pred1$pvals$predicted.value[1:3]
pred2se - pred1$pvals$std.error[1:3]
A - pred1$sed[1:3,1:3]

## calculate standard errors using AD chol

# use AD cholesky for se:
P <- list()
P[[1]] <- obj2$C
P[[2]] <- crossprod.spam(as.spam(Dg1))
P[[3]] <- crossprod.spam(as.spam(Dg2))
P[[4]] <- crossprod.spam(as.spam(Dg3))

C <- Reduce(`+`, P)
display(C)

ADobjC <- LMMsolver:::ADchol(P)

dim(C)
dlogdet1 <- LMMsolver:::dlogdet(ADobjC,c(1.0,0.0,0.0,0.0))
all.equal(pred2se, sqrt(dlogdet1[-1]))
sqrt(dlogdet1[2])

# use AD cholesky for sed:
P <- list()
P[[1]] <- obj2$C
P[[2]] <- crossprod.spam(as.spam(Dg12))
P[[3]] <- crossprod.spam(as.spam(Dg13))
P[[4]] <- crossprod.spam(as.spam(Dg23))

C <- Reduce(`+`, P)
display(C)

ADobjC <- LMMsolver:::ADchol(P)

dim(C)
dlogdet1 <- LMMsolver:::dlogdet(ADobjC,c(1.0,0.0,0.0,0.0))
A2 <- matrix(data=NA,ncol=3,nrow=3)
A2[1,2] <- A2[2,1] <- sqrt(dlogdet1[2])
A2[1,3] <- A2[3,1] <- sqrt(dlogdet1[3])
A2[2,3] <- A2[3,2] <- sqrt(dlogdet1[4])
all.equal(A, A2)

