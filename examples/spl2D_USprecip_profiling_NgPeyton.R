rm(list=ls())
#library(SAP)
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
grid <- c(300, 200)
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

dx1 <- attr(knots[[1]], which = 'dx')
dx2 <- attr(knots[[2]], which = 'dx')

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

C <- LMMsolver:::linearSum(theta, listP)
cholC <- chol(C)
slotNames(cholC)
sn <- cholC@supernodes
length(sn)
ns <- diff(sn, diff=1)
table(ns)

obj0 <- LMMsolver:::ADchol(listP)

det0 <- LMMsolver:::logdet(obj0, lambda=theta)
det1 <- 2.0*as.numeric(determinant(cholC)$modulus)

det0
det1

det0-det1

ED <- theta * LMMsolver:::dlogdet(obj0,theta)
ED

U <- cbind(obj1$X, obj1$Z)
tU <- t(U)
funSpam <- function(A, B) { A %*% B}

funMP <- function(A,B)
{
  LMMsolver:::MatrixProduct(A,B)
}

cntMP <- function(A,B)
{
  LMMsolver:::cntProduct(A,B)
}

cntMP(tU,U)
x1 <- funSpam(tU,U)
x2 <- funMP(tU,U)
all.equal(x1, x2)

# compare computation time
microbenchmark(funSpam(tU,U), funMP(tU,U), cntMP(tU, U), times=25L)


# compare direct way using dlogdet, and second option
#dlogdet1 <- as.numeric(LMMsolver:::dlogdet(obj0,theta))
#A <- LMMsolver:::DerivCholesky(cholC);
#dlogdet2 <- sapply(listP, function(x) {sum(A * x)})

#dlogdet1
#dlogdet2
#dlogdet1-dlogdet2

# compare computation times:

funlogdet <- function(lambda) {LMMsolver:::logdet(obj0, lambda)}
funSpam <- function(lambda) {
  C <- LMMsolver:::linearSum(lambda, listP)
  cholC <- update(cholC, C)
  2.0*as.numeric(determinant(cholC)$modulus)
}
funED <- function(lambda)
{
  ED <- lambda * LMMsolver:::dlogdet(ADobj,lambda)
  ED
}

# compare computation time
microbenchmark(funSpam(theta), funlogdet(theta), funED(theta), times=100L)

s <- proc.time()[3]
pred <- obtainSmoothTrend(obj1, grid=grid, includeIntercept = TRUE)
e <- proc.time()[3]

cat("Time pred: ", e-s, " seconds\n")

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

ggplot(plotDat, aes(x = lon, y = lat, fill = se)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradientn(colours = topo.colors(100))+
  labs(title = "Precipitation anomaly standard errors", x = "Longitude", y = "Latitude") +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

h1 <- (max(dat$lon) - min(dat$lon)) / nseg[1]
h2 <- (max(dat$lat) - min(dat$lat)) / nseg[2]
lambda1 <- as.numeric(obj1$sigma2e / obj1$tau2e[1]/ h1 ^ 3)
lambda2 <- as.numeric(obj1$sigma2e / obj1$tau2e[2]/ h2 ^ 3)
lambda1
lambda2





