## see also examples in ?gam.models (e.g. 'by' variables,
## random effects and tricks for large binary datasets)
rm(list=ls())
library(spam)
library(LMMsolver)
library(mgcv)

set.seed(2) ## simulate some data...
dat <- gamSim(1, n=400,dist="normal", scale=2)
nknots = 20

xmin = 0
xmax = 1
range = c(xmin,xmax)
b <- gam(y~s(x0,bs="ps",k=nknots) + s(x1,bs="ps",k=nknots)
         + s(x2,bs="ps",k=nknots) + s(x3,bs="ps",k=nknots),
         data=dat,method='REML',
         knots = list(x0 = range, x1=range, x2=range, x3=range))
summary(b)
plot(b,page=1)

nseg <- c(20, 50, 100,20)
# (solve in library)
eps <- 0.0001
knots0 <- PsplinesKnots(min(dat$x0)-eps,max(dat$x0)+eps,degree=3,nseg[1])
knots1 <- PsplinesKnots(min(dat$x1)-eps,max(dat$x1)+eps,degree=3,nseg[2])
knots2 <- PsplinesKnots(min(dat$x2)-eps,max(dat$x2)+eps,degree=3,nseg[3])
knots3 <- PsplinesKnots(min(dat$x3)-eps,max(dat$x3)+eps,degree=3,nseg[4])

B0 <- Bsplines(knots0, dat$x0)
B1 <- Bsplines(knots1, dat$x1)
B2 <- Bsplines(knots2, dat$x2)
B3 <- Bsplines(knots3, dat$x3)
q0 <- ncol(B0)
q1 <- ncol(B1)
q2 <- ncol(B2)
q3 <- ncol(B3)

Usc0 <- calcUsc(q0, 2)
Usc1 <- calcUsc(q1, 2)
Usc2 <- calcUsc(q2, 2)
Usc3 <- calcUsc(q3, 2)

lZ <- list()
lZ[[1]] = B0 %*% Usc0
lZ[[2]] = B1 %*% Usc1
lZ[[3]] = B2 %*% Usc2
lZ[[4]] = B3 %*% Usc3

Z <- as.matrix(do.call("cbind", lZ))
dim(Z)
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("B0","B1","B2","B3"))

obj1 <- LMMsolve(y~x0+x1+x2+x3, randomMatrices=lM, data=dat_ext,
                display=TRUE,monitor=TRUE)
round(obj1$ED,2)
obj1$logL

lZ <- list()
lZ[[1]] = B0
lZ[[2]] = B1
lZ[[3]] = B2
lZ[[4]] = B3

Z <- as.matrix(do.call("cbind", lZ))
dim(Z)
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("B0","B1","B2","B3"))

D0 <- diff(diag(q0),diff=2)
DtD0 <- crossprod(D0)
C0 <- matrix(data=0,ncol=2,nrow=q0)
C0[1,1] = C0[q0,2] = 1
P0 <- as.spam(DtD0 + tcrossprod(C0))

D1 <- diff(diag(q1),diff=2)
DtD1 <- crossprod(D1)
C1 <- matrix(data=0,ncol=2,nrow=q1)
C1[1,1] = C1[q1,2] = 1
P1 <- as.spam(DtD1 + tcrossprod(C1))

D2 <- diff(diag(q2),diff=2)
DtD2 <- crossprod(D2)
C2 <- matrix(data=0,ncol=2,nrow=q2)
C2[1,1] = C2[q2,2] = 1
P2 <- as.spam(DtD2 + tcrossprod(C2))

D3 <- diff(diag(q3),diff=2)
DtD3 <- crossprod(D3)
C3 <- matrix(data=0,ncol=2,nrow=q3)
C3[1,1] = C3[q3,2] = 1
P3 <- as.spam(DtD3 + tcrossprod(C3))

lGinv <- list()
lGinv[['B0']] <- P0
lGinv[['B1']] <- P1
lGinv[['B2']] <- P2
lGinv[['B3']] <- P3
names(lGinv)
obj2 <- LMMsolve(y~x0+x1+x2+x3, randomMatrices=lM,lGinverse=lGinv, data=dat_ext,
                display=TRUE,monitor=TRUE)
round(obj1$ED,2)
round(obj2$ED,2)
obj2$logL
obj1$logL-obj2$logL

# 2x4 = 8 sum to zero constraints:
cf <- coef(obj2)
cf$B0[c(1, q0)]
cf$B1[c(1, q1)]
cf$B2[c(1, q2)]
cf$B3[c(1, q3)]

