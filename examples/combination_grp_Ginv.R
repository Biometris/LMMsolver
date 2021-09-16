# Martin Boer, Biometris, WUR
#
# simple example of $a$ genotypes $n$ times replicated to illustrate
# how for defined matrices in the dataframe we can define a precision matrix
#
rm(list=ls())
library(asreml)
library(LMMsolver)
library(dplyr)
library(spam)

# a genotypes, n times replicated
a <- 10
n <- 3

set.seed(1234)

sigma2g <- 1.0
sigma2e <- 1.0

geno_eff  =rnorm(a,sd=sqrt(sigma2g))
geno_name = paste0("g",formatC(1:a,width=2,flag="0"))

dat <- data.frame(geno=rep(geno_name, each=n),
                  rep=rep(paste0("rep",1:n),times=a),
                  y = rep(geno_eff,each=n) + rnorm(a*n,sd=sigma2e),
                  stringsAsFactors = TRUE)

# LMMsolver:
obj1 <- LMMsolve(y~rep, random=~geno, data=dat, tolerance=1.0e-12,display=TRUE)
obj1$logL

# effective dimension:
obj1$ED

# alternative formulation, using transformation vector U = [1 D'],
#  y = Xa + Z v + e -> v = U w -> w = U^{-1} v, w ~ N(0,(UtU)^{-1}))
#
D <- diff(diag(a), diff=1)
U <- cbind(1, t(D))
UtU <- crossprod(U)

# here we use combination of argument group with
# lGinverse:
lZ <- list()
lZ[[1]] <- kronecker(U, rep(1,n))
Z <- do.call("cbind",lZ)
dat_ext = cbind(dat,Z)
lM <- ndxMatrix(dat, lZ, c("diff_geno"))
lGinv <- list()
lGinv[['diff_geno']] <- as.spam(UtU)

obj2 <- LMMsolve(y~rep, group=lM, lGinverse=lGinv,
                 data=dat_ext, tolerance=1.0e-12,display=TRUE)
obj1$logL
obj2$logL

obj1$ED
obj2$ED

coef(obj1)$'(Intercept)'
coef(obj2)$'(Intercept)'

coef(obj1)$rep
coef(obj2)$rep

a1 <- coef(obj1)$geno
a2 <- as.vector(U %*% coef(obj2)$diff_geno)
all.equal(a1, a2)

