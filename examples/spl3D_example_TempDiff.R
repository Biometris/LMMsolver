rm(list=ls())
library(mvtnorm)
library(SAP)
library(LMMsolver)

# testdata Daniela NPEC, see list of closed issues LMMsolver:
# "error after a number of iterations with spl3Dfast"
df <- read.csv("./examples/exampleTempDiff.csv")
head(df)

y <- df$par
x1 <- df$row
x2 <- df$col
x3 <- df$time

# set parameters:
knots <- c(2, 2, 50)
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

s <- proc.time()[3]
obj1 <- LMMsolve(fixed = par~1,
                 spline = ~spl3D(row, col, time, nseg=knots,
                              x1lim = x1lim, x2lim = x2lim, x3lim = x3lim),
                 data = df,
                 trace = trace,
                 tolerance = thr)
e <- proc.time()[3]
cat("Proc time", e-s, " seconds\n")
summary(obj1)

# compare effective dimensions:
obj0$edf  # SAP package
obj1$ED   # LMMsolve with spatial argument

# compare smooth trends
fit1 <- obtainSmoothTrend(obj1, grid=grid, includeIntercept = TRUE)
range(fit0 - fit1$ypred)

## Obtain smooth trend for row 0, col 0 only.
newdat <- df[df$row == 0 & df$col == 0, ]
fit2 <- obtainSmoothTrend(obj1, newdata = newdat, includeIntercept = TRUE)


