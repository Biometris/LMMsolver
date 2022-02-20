# Simulate data, using first example JOPS book:
library(ggplot2)
library(LMMsolver)
library(spam)
library(agridat)

n = 150
set.seed(2016)
x = seq(0, 1, length = n)
sigma2e = 0.15^2
simFun <- function(x) {0.3 + sin(1.2 * x + 0.3) }

y =  simFun(x) + rnorm(n) * sqrt(sigma2e)
dat <- data.frame(x=x,y=y)

# Make a matrix containing the B-spline basis
nseg = 15

# solve
obj <-LMMsolve(y~1, spline=~spl1D(x,nseg=nseg, scaleX=FALSE, xlim=c(0,1),
                                  pord=2),data=dat)
summary(obj)

nGrid <- 1000
plotDat <- obtainSmoothTrend(obj, grid=nGrid, includeIntercept = TRUE)
x0 <- seq(0, 1, length = nGrid)
plotDat$ysim <- simFun(x0)

ggplot(data = dat, aes(x = x, y = y)) +
  geom_point(size = 1.2) +
  geom_line(data = plotDat, aes(y = ypred), color = "red", size = 1) +
  geom_line(data = plotDat, aes(y = ysim), color = "green", size = 1) +
  geom_line(data = plotDat, aes(y = ypred-2*se), col='blue', size=1) +
  geom_line(data = plotDat, aes(y = ypred+2*se), col='blue', size=1) +
  theme(panel.grid = element_blank())

# use eigendecomposition
knots <- LMMsolver:::PsplinesKnots(xmin=0,xmax=1,degree=3,nseg=nseg)
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
