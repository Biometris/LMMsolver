library(LMMsolver)
library(spam)

# Get precipitation data from spam
data(USprecip)
dat = data.frame(USprecip)

datOrig <- dat

# only use observed data
dat = subset(dat, infill==1)
nrow(dat) # 5906 true records, as in SAP2014 paper.

# arbitrary split of US in west and east, to test conditional factor...
dat$Region <- as.factor(ifelse(dat$lon < -100, "west", "east"))

s <- proc.time()[3]
obj1 <- LMMsolve(fixed = anomaly~Region,
                 spline = ~spl2D(x1 = lon, x2 = lat, nseg = c(20, 20),
                                  , cond=Region, level="west") +
                           spl2D(x1 = lon, x2 = lat, nseg = c(20, 20),
                                    cond=Region, level="east"),
                 data = dat,
                 trace = TRUE)
e <- proc.time()[3]
cat("Time LMMsolve ", e-s, " seconds\n")
summary(obj1)

