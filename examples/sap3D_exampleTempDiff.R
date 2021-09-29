rm(list=ls())
library(mvtnorm)
library(SAP)
library(LMMsolver)

# testdata Daniela NPEC, see list of closed issues LMMsolver:
# "error after a number of iterations with spl3Dfast"
df <- read.csv("examples/exampleTempDiff.csv")
head(df)

y <- df$par
x1 <- df$row
x2 <- df$col
x3 <- df$time

# new data for prediction:
#nNew <- 5
#newData <- data.frame(x1=runif(nNew),x2=runif(nNew),x3=runif(nNew))

#
# set parameters:
#
knots <- c(2, 2, 50)
grid <- c(30,40,50)
trace <- TRUE
thr <- 1.0e-8  # convergence tolerance
######################

# original sap package:
obj0 <- SAP::sap3D(y, x1, x2, x3, knots=knots, trace=trace, thr=thr)
obj0$edf

fit0 <- predict(obj0, grid=grid)$eta
#fit0n <- predict(obj0, newdata=newData)$eta

# fast sap, using LMMsolver:
obj1 <- spl3Dfast(y, x1, x2, x3, nseg=knots, trace=trace, tolerance=thr)
obj1$edf
fit1 <- predict(obj1, grid=grid)$eta

# extra space, to be consitent with knot positions in SAP package:
x1lim <- c(min(x1)-0.01, max(x1)+0.01)
x2lim <- c(min(x2)-0.01, max(x2)+0.01)
x3lim <- c(min(x3)-0.01, max(x3)+0.01)

obj2 <- LMMsolve(fixed = par~1,
                 spline = ~LMMsolver::spl3D(row, col, time, nseg=knots,
                              x1lim = x1lim, x2lim = x2lim, x3lim = x3lim),
                 data = df,
                 trace = trace,
                 tolerance = thr)
summary(obj2)

#fit1n <- predict(obj1, newdata=newData)$eta

# fast sap, using LMMsolver:
#obj2 <- spl3Dfast(y, x1, x2, x3, knots=knots, trace=trace, thr=thr,scaleX=TRUE)
#obj2$edf

#fit2 <- predict(obj2, grid=grid)$eta
#fit2n <- predict(obj2, newdata=newData)$eta

# compare fit newdata:
#range(fit1n-fit0n)

# compare effective dimensions:
obj0$edf  # SAP package
obj1$edf  # spl3D fast function
obj2$ED   # LMMsolve with spatial argument

fit2 <- obtainSmoothTrend3D(obj2, grid=grid)$eta

# compare fit on grid:
range(fit1-fit0)
range(fit2-fit0)
range(fit1-fit2)
