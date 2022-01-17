rm(list=ls())
library(SAP)
library(LMMsolver)
library(spam)
library(fields)
library(ggplot2)
library(maps)
library(dplyr)
library(microbenchmark)

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

x1lim <- c(min(x1)-0.01, max(x1)+0.01)
x2lim <- c(min(x2)-0.01, max(x2)+0.01)

#
# use spatial option in LMMsolve:
#
s <- proc.time()[3]
obj1 <- LMMsolve(fixed = anomaly~1,
                 spline = ~spl2D(x1 = lon, x2 = lat, nseg = nseg,
                                  x1lim=x1lim, x2lim=x2lim),
                 data = dat,
                 trace = trace,
                 tolerance = thr)
e <- proc.time()[3]
cat("Time LMMsolve ", e-s, " seconds\n")
EDtbl <- summary(obj1)

# compare effective dimensions, in SAP2014 paper
# ED(lon) = 302.656
# ED(lat) = 408.757
EDtbl

knots <- list()
knots[[1]] <- LMMsolver:::PsplinesKnots(x1lim[1], x1lim[2], degree = 3, nseg = nseg[1])
knots[[2]] <- LMMsolver:::PsplinesKnots(x2lim[1], x2lim[2], degree = 3, nseg = nseg[2])

x1 = dat$lon
x2 = dat$lat
B1 <- LMMsolver:::Bsplines(knots[[1]], x1)
B2 <- LMMsolver:::Bsplines(knots[[2]], x2)
q <- c(ncol(B1), ncol(B2))

B <- LMMsolver:::RowKronecker(B1, B2)

P1 <- LMMsolver:::constructPenalty(q[1],2) %x% diag(1,q[2])
P2 <- diag(1,q[1]) %x% LMMsolver:::constructPenalty(q[2],2)

theta <- EDtbl$Penalty[c(5,3,4)]
listP <- list()
listP[[1]] <- as.spam(crossprod(B))
listP[[2]] <- as.spam(P1)
listP[[3]] <- as.spam(P2)

sapply(listP,dim)

ADobj <- LMMsolver:::ADchol(listP)

EDtbl

ED <- theta * LMMsolver:::dlogdet(ADobj,theta)

# supernodes....
displayMME(obj1, cholesky=TRUE)

C <- theta[1]*listP[[1]] + theta[2]*listP[[2]] + theta[3]*listP[[3]]
cholC <- chol(C)
slotNames(cholC)
sn <- cholC@supernodes
length(sn)
ns <- diff(sn, diff=1)
table(ns)

obj0 <- LMMsolver:::ADchol(listP)
obj1 <- LMMsolver:::ADcholNgPeyton(listP)

det0 <- LMMsolver:::logdet(obj0, lambda=theta)
det1 <- LMMsolver:::logdetNgPeyton(obj1, lambda=theta)
det2 <- 2.0*as.numeric(determinant(cholC)$modulus)

ED2 <- theta * LMMsolver:::dlogdetNgPeyton(obj1,theta)
ED
ED2

det0
det1
det2

det0-det1
det0-det2
det1-det2

funOrg <- function(lambda) {LMMsolver:::logdet(obj0, lambda)}
funNew <- function(lambda) {LMMsolver:::logdetNgPeyton(obj1, lambda)}
funSpam <- function(lambda) {
  C=lambda[1]*listP[[1]] + lambda[2]*listP[[2]] + lambda[3]*listP[[3]]
  cholC <- update(cholC, C)
  2.0*as.numeric(determinant(cholC)$modulus)
}
funED <- function(lambda)
{
  ED <- lambda * LMMsolver:::dlogdet(ADobj,lambda)
  ED
}

funED2 <- function(lambda)
{
  ED <- lambda * LMMsolver:::dlogdetNgPeyton(obj1,lambda)
  ED
}

# compare computation time
microbenchmark(funSpam(theta), funOrg(theta), funNew(theta), funED(theta),
            funED2(theta), times=10L)



