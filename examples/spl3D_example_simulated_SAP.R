rm(list=ls())
library(mvtnorm)
library(SAP)
library(LMMsolver)

set.seed(1234)

Wood3D <- function(x1, x2, x3) { # Wood 2006
  y = 1.5*exp(-(x1-0.2)^2/5 -(x2-0.5)^2/3 -(x3-0.9)^2/4)
    + 0.5*exp(-(x1-0.3)^2/4 -(x2-0.7)^2/2 - (x3-0.4)^2/6)
    + exp(-(x1-0.1)^2/5 -(x2-0.3)^2/5 - (x3-0.7)^2/4)
  y
}

N <- 500
x1 <- runif(N)
x2 <- runif(N)
x3 <- runif(N)
eps <- rnorm(N,sd=0.1)

y = Wood3D(x1, x2, x3) + eps

# set parameters:
knots <- c(4, 4, 4)
grid <- c(30,40,50)
trace <- TRUE
thr <- 1.0e-8  # convergence tolerance
######################

# original sap package:
obj0 <- SAP::sap3D(y, x1, x2, x3, knots=knots, trace=trace, thr=thr)
obj0$edf

fit0 <- predict(obj0, grid=grid)$eta

# extra space, to be consistent with knot positions in SAP package:
x1lim <- c(min(x1)-0.01, max(x1)+0.01)
x2lim <- c(min(x2)-0.01, max(x2)+0.01)
x3lim <- c(min(x3)-0.01, max(x3)+0.01)

df <- data.frame(y=y,x1=x1,x2=x2,x3=x3)
obj1 <- LMMsolve(fixed = y~1,
                 spline = ~spl3D(x1, x2, x3, nseg=knots,
                              x1lim = x1lim, x2lim = x2lim, x3lim = x3lim),
                 data = df,
                 trace = trace,
                 tolerance = thr)
summary(obj1)

# compare effective dimensions:
obj0$edf  # SAP package
obj1$ED   # LMMsolve with spatial argument

# compare smooth trends
fit1 <- obtainSmoothTrend(obj1, grid=grid, includeIntercept = TRUE)
range(fit0 - fit1$ypred)

