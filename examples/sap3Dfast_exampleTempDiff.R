rm(list=ls())
library(mvtnorm)
library(SAP)
library(LMMsolver)


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
obj0 <- sap3D(y, x1, x2, x3, knots=knots, trace=trace, thr=thr)
obj0$edf

fit0 <- predict(obj0, grid=grid)$eta
#fit0n <- predict(obj0, newdata=newData)$eta

# fast sap, using LMMsolver:
obj1 <- sap3Dfast(y, x1, x2, x3, knots=knots, trace=trace, tolerance=thr)
obj1$edf
fit1 <- predict(obj1, grid=grid)$eta

obj2 <- LMMsolve(fixed = formula(y~1),
                 spatial = ~sap3D(x1, x2, x3, knots),
                 data = df,
                 trace = trace,
                 tolerance = thr)
obj2$ED


#fit1n <- predict(obj1, newdata=newData)$eta

# fast sap, using LMMsolver:
#obj2 <- sap3Dfast(y, x1, x2, x3, knots=knots, trace=trace, thr=thr,scaleX=TRUE)
#obj2$edf

#fit2 <- predict(obj2, grid=grid)$eta
#fit2n <- predict(obj2, newdata=newData)$eta

# compare fit on grid:
range(fit1-fit0)

# compare fit newdata:
#range(fit1n-fit0n)

# compare effective dimensions:
obj0$edf - obj1$edf
