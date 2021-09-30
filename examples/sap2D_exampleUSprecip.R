rm(list=ls())
library(SAP)
library(LMMsolver)
library(spam)
library(fields)

# Get precipitation data from spam
data(USprecip)
dat = data.frame(USprecip)
# only use observed data
dat = subset(dat, infill==1)
nrow(dat) # 5906 true records, as in SAP2014 paper.

y <- dat$anomaly
x1 <- dat$lon
x2 <- dat$lat

#
# set parameters:
#
innerKnots <- 40 # See SAP 2014 paper
knots <- nseg <- c(innerKnots+1, innerKnots+1)
grid <- c(300,200)
trace <- TRUE
thr <- 1.0e-7  # convergence tolerance
######################

# original sap package:
obj0 <- SAP::sap2D(y, x1, x2 , knots=knots, trace=trace, thr=thr)
obj0$edf
fit0 <- predict(obj0, grid=grid)$eta
# reorder fit....
fit0 <- matrix(data=fit0, nrow=grid[1], ncol=grid[2], byrow=FALSE)

# fast sap, using LMMsolver:
obj1 <- spl2Dfast(y, x1, x2, nseg=knots, trace=trace, tolerance=thr,scaleX=FALSE)
obj1$edf
fit1 <- predict(obj1, grid=grid)$eta
fit1 <- matrix(data=fit1, nrow=grid[1], ncol=grid[2], byrow=TRUE)

# fast sap, using LMMsolver:
obj2 <- spl2Dfast(y, x1, x2, nseg=knots, trace=trace, tolerance=thr,scaleX=TRUE)
obj2$edf
fit2 <- predict(obj2, grid=grid)$eta
fit2 <- matrix(data=fit2, nrow=grid[1], ncol=grid[2], byrow=TRUE)

# extra space, to be consistent with knot positions in SAP package:
x1lim <- c(min(x1)-0.01, max(x1)+0.01)
x2lim <- c(min(x2)-0.01, max(x2)+0.01)

#
# use spatial option in LMMsolve:
#
obj3 <- LMMsolve(fixed = anomaly~1,
                 spline = ~LMMsolver::spl2D(x1 = lon, x2 = lat, nseg = nseg,
                              scaleX = FALSE, x1lim=x1lim, x2lim=x2lim),
                 data = dat,
                 trace = trace,
                 tolerance = thr)

obj4 <- LMMsolve(fixed = anomaly~1,
                 spline = ~LMMsolver::spl2D(x1 = lon, x2 = lat, nseg = nseg,
                                             scaleX = TRUE,
                                            x1lim=x1lim, x2lim=x2lim),
                 data = dat,
                 trace = trace,
                 tolerance = thr)
summary(obj4)

# compare effective dimensions, in SAP2014 paper
# ED(lon) = 302.656
# ED(lat) = 408.757
obj0$edf
obj1$edf
obj2$edf
obj3$ED
obj4$ED

fit3 <- LMMsolver:::obtainSmoothTrend2D(obj3, grid)$eta
fit3a <- obtainSmoothTrend(obj3, grid)$eta

fit3 <- matrix(data=fit3, nrow=grid[1], ncol=grid[2], byrow=TRUE)

fit4 <- obtainSmoothTrend(obj4, grid)$eta
fit4 <- matrix(data=fit4, nrow=grid[1], ncol=grid[2], byrow=TRUE)

# compare fit on grid:
range(fit0 - fit1)
range(fit0 - fit2)
range(fit0 - fit3)
range(fit0 - fit4)

