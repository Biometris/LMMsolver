library(SOP)
library(spam)
library(LMMsolver)
library(ggplot2)
library(JOPS)
library(gridExtra)
library(colorspace)

savePlotTimeTable <- TRUE

# Get precipitation data from spam
data(USprecip)
dat <- data.frame(USprecip)
# only use observed data
dat <- subset(dat, infill==1)
nrow(dat) # 5906 true records, as in SAP2014 paper.

y <- dat$anomaly
x1 <- dat$lon
x2 <- dat$lat

nrSegments <- seq(20, 80, by=5)
K <- length(nrSegments)
LMMsolver_time <- rep(NA,K)
#SOP_time <- rep(NA,K)
#ED_SOP <- matrix(data=NA,ncol=2,nrow=K)
ED_LMMsolver <- matrix(data=NA,ncol=2,nrow=K)
s <- proc.time()[3]
for (i in 1:K) {
  curTime <- proc.time()[3]
  cat("Nseg", nrSegments[i], " time", curTime - s, "\n")
  nseg <- c(nrSegments[i], nrSegments[i])

  s1 <- proc.time()[3]
  obj1 <- LMMsolve(fixed = anomaly~1,
                 spline = ~spl2D(x1 = lon, x2 = lat, nseg = nseg),
                 data = dat)
  e1 <- proc.time()[3]
  LMMsolver_time[i] <- e1-s1
}
LMMsolver_time

df_times <- data.frame(nseg=nrSegments, time = LMMsolver_time)
df_times
