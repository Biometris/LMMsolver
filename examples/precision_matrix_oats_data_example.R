# Martin Boer, Biometris, Wageningen

# example of definition precision matrix
# with LMMsolver
#
rm(list=ls())

library(agridat)
library(LMMsolver)
library(spam)
library(splines)
library(dplyr)
library(asreml)

# oats data
n <- 4  # number of units per block
b <- 6  # number of blocks per Rep
r <- 3  # number of Reps
v <- 24 # number of genotypes/Rep:
N <- n*b*r # total number of observations.

data(john.alpha)
dat <- john.alpha
dat$plot_in_rep = as.factor((dat$plot - 1) %% v + 1)
dat$plot = as.factor(dat$plot)

dat <- rename(dat, Trt = gen)
dat <- rename(dat, Block = block)
dat <- rename(dat, Rep = rep)
dat <- rename(dat, Row=plot)
dat <- rename(dat, rRow = plot_in_rep)
dat <- rename(dat, y=yield)
oats.df <- dat
rm(dat)

# define a precision matrix for Linear variance model, or Random walk:
invM_matrix <- function(n, type) {
  D = diff(diag(n), diff=1)
  Delta = 0.5*crossprod(D)
  if (type=='RWf') {
    alpha <- 1.0
    Con <- c(1,rep(0,n-1))
  } else if (type=='RWb') {
    alpha <- 1.0
    Con <- c(rep(0,n-1),1)
  } else if (type=="LV") {
    alpha <- 0.5
    Con <- c(sqrt(1/(n-1)),rep(0,n-2),sqrt(1/(n-1)))
  } else
  {
    stop("invM type not not defined")
  }
  M <- alpha*(2*Delta + tcrossprod(Con))
  rownames(M) <- paste0(1:n)
  colnames(M) <- colnames(1:n)
  attr(M, 'INVERSE') <- TRUE
  M
}

# test with asreml using linear variance model:
Ginv_LV2 <- invM_matrix(v, "LV")
asr = asreml(y ~ Trt + Rep, random = ~idv(Rep):vm(rRow, Ginv_LV2),
              data = oats.df,trace=FALSE, maxiter=25)

# compare with LMMsolve using inverse matrix:
Ginv <- as.spam(Ginv_LV2)
lGinv <- list()
lGinv[["rRow:Rep"]] = kronecker.spam(diag(r), Ginv)
obj1 <- LMMsolve(y~Trt+Rep, random=~rRow:Rep, lGinverse=lGinv,data=oats.df, eps=1.0e-10)

# check the buildin-constraints:
LV_eff <- coef(obj1)$'rRow:Rep'
LV_eff[c(1,v+1,2*v+1)] + LV_eff[c(v,2*v,3*v)]

# P-splines formulation:
lZ <- list()
Usc <- calcUsc(v, 1)
lZ[[1]] = kronecker(diag(1,r), Usc)
Z <- do.call("cbind",lZ)
oats.df_ext = cbind(oats.df,Z)
lM <- ndxMatrix(oats.df, lZ, c("P-splines"))
obj2 <- LMMsolve(y~Trt+Rep, group=lM, data=oats.df_ext,monitor=FALSE,eps=1.0e-10)

D = diff(diag(v), diff=1)
DtD = crossprod(D)
P = kronecker.spam(diag(r), DtD)
d <- eigen(P)$values
U0 <- eigen(P)$vectors[,which(abs(d)<1.0e-10)]
dim(U0)

# Use as constraint U0
C <- U0
CtG <- t(C) %*% U0
CtG
lGinv <- list()
lGinv[["rRow:Rep"]] = as.spam(P + tcrossprod(C))
obj3 <- LMMsolve(y~Trt+Rep, random=~rRow:Rep, lGinverse=lGinv,
                 data=oats.df, eps=1.0e-10, display=TRUE)
p <- v + r -1
random_eff <- obj3$a[-c(1:p)]

# built-in constraints:
t(U0) %*% random_eff

# automatic selection of constraint C based on null space:

# ndx for sparse constraint, it doesn't matter
# which one we choose...
set.seed(1234)
ndxCon <- rep(NA, r)
C <- matrix(data=0,nrow=N,ncol=3)
for (i in 1:r)
{
  ndx_non_zero <- which(U0[,i]!=0)
  # choose the next element, not overlapping with previous indices:
  k <- sample(setdiff(ndx_non_zero, ndxCon), 1)
  C[k,i] = 1/U0[k,i]
  ndxCon[i] = k
}
ndxCon

# CtG is identity matrix
CtG <- t(C) %*% U0
CtG
lGinv <- list()
lGinv[["rRow:Rep"]] = as.spam(P + tcrossprod(C))
obj4 <- LMMsolve(y~Trt+Rep, random=~rRow:Rep, lGinverse=lGinv,
                 data=oats.df, eps=1.0e-10, display=TRUE)
p <- v + r -1
random_eff <- obj4$a[-c(1:p)]
random_eff[ndxCon]

# compare the results, should be identical:

logL <- rep(NA,5)
logL[1] = asr$loglik
logL[2] = obj1$logL
logL[3] = obj2$logL
logL[4] = obj3$logL
logL[5] = obj4$logL

# compare all deviances:
abs(apply(combn(logL,2), 2, diff))

obj1$ED
obj2$ED
obj3$ED
obj4$ED

cf1 <- coef(obj1)
cf2 <- coef(obj2)
cf3 <- coef(obj3)
cf4 <- coef(obj4)

cf1$'(Intercept)'
cf2$'(Intercept)'
cf3$'(Intercept)'
cf4$'(Intercept)'

cf1$Rep
cf2$Rep
cf3$Rep
cf4$Rep

range(cf1$Trt - cf2$Trt)
range(cf1$Trt - cf3$Trt)
range(cf1$Trt - cf4$Trt)


