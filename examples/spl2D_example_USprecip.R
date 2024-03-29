rm(list=ls())
#library(SAP)
library(LMMsolver)
library(spam)
library(fields)
library(ggplot2)
library(maps)
library(dplyr)
library(SOP)

# Get precipitation data from spam
data(USprecip)
dat = data.frame(USprecip)

datOrig <- dat

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

# sop package
obj0 <- sop(formula = anomaly ~ f(x1=lon, x2=lat, nseg = nseg), data = dat,
          control = list(trace = TRUE, epsilon=thr))
obj0$deviance

# original sap package:
#obj0 <- SAP::sap2D(y, x1, x2 , knots=knots, trace=trace, thr=thr)
#obj0$edf

obj0$f$`f(x1 = lon, x2 = lat, nseg = nseg)`$Xmat$terms$MM$MM1$knots
range(x1)

obj0$f$`f(x1 = lon, x2 = lat, nseg = nseg)`$Xmat$terms$MM$MM2$knots
range(x2)

newdat <- expand.grid(lon=seq(min(x1),max(x1),length=grid[1]),
                       lat=seq(min(x2),max(x2),length=grid[2]))

fit0 <- predict(obj0, newdata=newdat)
# reorder fit....
Mfit0 <- matrix(data=fit0, nrow=grid[1], ncol=grid[2], byrow=FALSE)

# extra space, to be consistent with knot positions in SAP package:
x1lim <- c(min(x1), max(x1))
x2lim <- c(min(x2), max(x2))

#
# use spatial option in LMMsolve:
#
s <- proc.time()[3]
obj1 <- LMMsolve(fixed = anomaly~1,
                 spline = ~spl2D(x1 = lon, x2 = lat, nseg = nseg),
                #                  x1lim=x1lim, x2lim=x2lim),
                 data = dat,
                 trace = trace,
                 tolerance = thr)
e <- proc.time()[3]
cat("Time LMMsolve ", e-s, " seconds\n")
summary(obj1)

# compare effective dimensions, in SAP2014 paper
# ED(lon) = 302.656
# ED(lat) = 408.757
obj0$out$edf
obj1$ED

# only use this in stead of grid to keep consistent with SAP
#grid1 <- seq(min(x1), max(x1),  length=grid[1])
#grid2 <- seq(min(x2), max(x2),  length=grid[2])
#newdata = expand.grid(grid1, grid2)
#head(newdata)
#colnames(newdata) = c('lon','lat')

pred <- obtainSmoothTrend(obj1, newdata=newdat, includeIntercept = TRUE)
Mfit1 <- matrix(data=pred$ypred, nrow=grid[1], ncol=grid[2], byrow=FALSE)

# compare fit on grid:
range(Mfit0 - Mfit1)

# make a plot
plotDat <- pred

usa = maps::map("usa", regions = "main", plot = FALSE)

v <- sp::point.in.polygon(plotDat$lon, plotDat$lat, usa$x, usa$y)

plotDat <- plotDat[v == 1, ]

ggplot(plotDat, aes(x = lon, y = lat, fill = ypred)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradientn(colours = topo.colors(100))+
  labs(title = "Precipitation anomaly", x = "Longitude", y = "Latitude") +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# predictions for new data, using coordinates
# from maps library:
newdat <- us.cities %>% dplyr::select(name, lat, long)
newdat <- newdat %>% rename(lon=long)
head(newdat)

pred2 <- obtainSmoothTrend(obj1, newdata=newdat, includeIntercept = TRUE)
head(pred2)

