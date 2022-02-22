# Simulated data use Jullion and Lambert 2007:

library(ggplot2)
library(LMMsolver)
library(spam)

set.seed(12345)

n <- 250
xmin <- -2
xmax <-  2
x <- seq(xmin, xmax, length = n)

# see Jullion and Lambert 2007:
simFun <- function(x) {
  (1+exp(-4*(x-0.3)))^(-1) +
  (1+exp( 3*(x-0.2)))^(-1) +
  (1+exp(-4*(x-0.7)))^(-1) +
  (1+exp( 5*(x-0.8)))^(-1)
}

sigma2e <- 0.0169

y <-  simFun(x) + rnorm(n, sd = sqrt(sigma2e))
dat <- data.frame(x=x,y=y)

# Make a matrix containing the B-spline basis
nseg <- 1000

# solve
obj <-LMMsolve(y~1, spline=~spl1D(x,nseg=nseg, scaleX=FALSE, xlim=c(xmin,xmax),
                                  pord=2),data=dat, trace=TRUE)
summary(obj)

nGrid <- 1000
plotDat <- obtainSmoothTrend(obj, grid=nGrid, includeIntercept = TRUE)
x0 <- seq(xmin, xmax, length = nGrid)
plotDat$ysim <- simFun(x0)

ggplot(data = dat, aes(x = x, y = y)) +
  geom_point(size = 1.2) +
  geom_line(data = plotDat, aes(y = ypred), color = "red", size = 1) +
  geom_line(data = plotDat, aes(y = ysim), color = "green", size = 1) +
  geom_line(data = plotDat, aes(y = ypred-2*se), col='blue', size=1) +
  geom_line(data = plotDat, aes(y = ypred+2*se), col='blue', size=1) +
  theme(panel.grid = element_blank())

# use eigendecomposition
knots <- LMMsolver:::PsplinesKnots(xmin=xmin,xmax=xmax,degree=3,nseg=nseg)
B <- LMMsolver:::Bsplines(knots, x)
q <- ncol(B)
dx <- attr(knots,which='dx')
DtDsc <- LMMsolver:::constructPenalty(q, pord=2, dx)
eig <- eigen(DtDsc)
Usc <- eig$vectors[, -c(q-1,q)] %*% diag(eig$values[-c(q-1,q)]^(-1/2))
Z <- B %*% Usc
dat_ext <- cbind(dat, Z)
L <- list(ndx = c(3:q))
obj2 <- LMMsolve(y~x,random=~grp(ndx), group=L, data=dat_ext)

obj$logL - obj2$logL

summary(obj)
summary(obj2)
