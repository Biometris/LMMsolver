rm(list=ls())
library(SAP)
library(LMMsolver)

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
knots <- c(41,41)
grid <- c(40,40)
trace <- TRUE
thr <- 1.0e-7  # convergence tolerance
######################

# original sap package:
obj0 <- sap2D(y, x1, x2 , knots=knots, trace=trace, thr=thr)
obj0$edf
fit0 <- predict(obj0, grid=grid)$eta
# reorder fit....
fit0 <- matrix(data=fit0, nrow=grid[1], ncol=grid[2], byrow=TRUE)

# fast sap, using LMMsolver:
obj1 <- sap2Dfast(y, x1, x2, knots=knots, trace=trace, tolerance=thr)
obj1$edf
fit1 <- predict(obj1, grid=grid)$eta

# fast sap, using LMMsolver:
obj2 <- sap2Dfast(y, x1, x2, knots=knots, trace=trace, tolerance=thr,scaleX=TRUE)
obj2$edf
fit2 <- predict(obj2, grid=grid)$eta

# compare fit on grid:
range(fit1-fit0)
range(fit2-fit0)
range(fit1-fit2)

# compare effective dimensions:
obj0$edf - obj1$edf
obj0$edf - obj2$edf
obj1$edf - obj2$edf

#
# use new spatial option in LMMsolve:
#
df = data.frame(y=y, x1=x1, x2=x2)

obj3 <- LMMsolve(fixed = formula(y~1),
                 spatial = ~LMMsolver::sap2D(x1, x2, knots),
                 data = df,
                 trace = trace,
                 tolerance = thr)

# effective dimensions for model 3 slightly different:
obj0$edf
obj1$edf
obj2$edf
obj3$ED

