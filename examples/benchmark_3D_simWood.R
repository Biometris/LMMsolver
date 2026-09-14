rm(list=ls())
library(SOP)
#library(spam)
library(LMMsolver)
library(ggplot2)

set.seed(1234)

Wood3D <- function(x1, x2, x3) { # Wood 2006
  y = 1.5*exp(-(x1-0.2)^2/5 -(x2-0.5)^2/3 -(x3-0.9)^2/4)
  + 0.5*exp(-(x1-0.3)^2/4 -(x2-0.7)^2/2 - (x3-0.4)^2/6)
  + exp(-(x1-0.1)^2/5 -(x2-0.3)^2/5 - (x3-0.7)^2/4)
  y
}

N <- 5000
x1 <- runif(N)
x2 <- runif(N)
x3 <- runif(N)
eps <- rnorm(N,sd=0.1)

y = Wood3D(x1, x2, x3) + eps

dat <- data.frame(y=y,x1=x1,x2=x2,x3=x3)

nrSegments <- seq(4, 12, by=1)
K <- length(nrSegments)
LMMsolver_time <- rep(NA,K)
ED_LMMsolver <- matrix(data=NA,ncol=3,nrow=K)
s <- proc.time()[3]
for (i in 1:K) {
  curTime <- proc.time()[3]
  cat("Nseg", nrSegments[i], " time", curTime - s, "\n")
  nseg <- rep(nrSegments[i],3)

  s1 <- proc.time()[3]
  obj1 <- LMMsolve(fixed = y~1,
                 spline = ~spl3D(x1 = x1, x2 = x2, x3=x3, nseg = nseg),
                 data = dat)
  e1 <- proc.time()[3]
  LMMsolver_time[i] <- e1-s1


}
packageVersion("LMMsolver")
df_times <- data.frame(nseg=nrSegments, time = LMMsolver_time)
df_times
