## see also examples in ?gam.models (e.g. 'by' variables,
## random effects and tricks for large binary datasets)
rm(list=ls())
library(spam)
library(LMMsolver)
library(mgcv)

set.seed(2) ## simulate some data...
dat <- gamSim(1, n=400, dist="normal", scale=2)
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

nseg <- c(20, 20, 80,20)

lKnots <- lapply(nseg, function(x) {PsplinesKnots(xmin, xmax, degree=3,x)})

# list of B-splines bases:
lB <- list()
lB[[1]] <- Bsplines(lKnots[[1]], dat$x0)
lB[[2]] <- Bsplines(lKnots[[2]], dat$x1)
lB[[3]] <- Bsplines(lKnots[[3]], dat$x2)
lB[[4]] <- Bsplines(lKnots[[4]], dat$x3)
q <- sapply(lB, ncol)
q

lUsc <- lapply(q, function(x) { calcUsc(x,2)})
lZ <- lapply(1:4,function(x) {lB[[x]] %*% lUsc[[x]]})
Z <- as.matrix(do.call("cbind", lZ))
dim(Z)
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("B0","B1","B2","B3"))

obj1 <- LMMsolve(y~x0+x1+x2+x3, group=lM, data=dat_ext, display=FALSE,monitor=TRUE)
round(obj1$ED,2)
obj1$logL

# sparse solution:
lZ <- lB
Z <- as.matrix(do.call("cbind", lZ))
dim(Z)
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("B0","B1","B2","B3"))

# make penalty matrix, with second order diff
# plus boundary constraints:
makeP <- function(q)
{
  D <- diff(diag(q), diff=2)
  DtD <- crossprod(D)
  C <- matrix(data=0, ncol=2, nrow=q)
  C[1,1] = C[q,2] = 1
  P <- as.spam(DtD + tcrossprod(C))
  P
}

lGinv <- lapply(q, makeP)
sapply(lGinv, dim)
names(lGinv) <- paste0("B",0:3)
names(lGinv)
obj2 <- LMMsolve(y~x0+x1+x2+x3, group=lM,lGinverse=lGinv, data=dat_ext,
                display=TRUE,monitor=TRUE)
round(obj1$ED,2)
round(obj2$ED,2)

obj2$logL
obj1$logL - obj2$logL

# 2x4 = 8 sum to zero constraints:
cf <- coef(obj2)
cf$B0[c(1, q[1])]
cf$B1[c(1, q[2])]
cf$B2[c(1, q[3])]
cf$B3[c(1, q[4])]

