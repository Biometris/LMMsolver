# Simulate data, using first example JOPS book:
library(ggplot2)
library(LMMsolver)

n = 150
set.seed(2016)
x = seq(0, 1, length = n)
sigma2e = 0.15^2

simFun <- function(x) {0.3 + sin(1.2 * x + 0.3) }

y =  simFun(x) + rnorm(n) * sqrt(sigma2e)
dat <- data.frame(x=x,y=y)

# Make a matrix containing the B-spline basis
nseg = 15

# solve
obj <-LMMsolve(y~1,spline=~spl1D(x,nseg=nseg, scaleX=FALSE, xlim=c(0,1)),data=dat)
summary(obj)

plotDat <- obtainSmoothTrend(obj, grid = 1000, includeIntercept = TRUE)

knots <- LMMsolver:::PsplinesKnots(xmin=0,xmax=1,degree=3,nseg=nseg)
B <- LMMsolver:::Bsplines(knots, x)
q <- ncol(B)
DtD <- LMMsolver:::constructPenalty(q, pord=2)
BtB <- crossprod(B)

sigma2e <- obj$sigma2e
lambda <- obj$theta[1]/obj$theta[2]

C <- (BtB+lambda*DtD)
BtY <- t(B) %*% y

lC <- list()
lC[[1]] <- BtB
lC[[2]] <- DtD
ADcholC <- LMMsolver:::ADchol(lC)

theta <- obj$theta[c(2,1)]

ED <- theta * LMMsolver:::dlogdet(ADcholC, theta)

sparseInverse <- LMMsolver:::DerivCholesky(chol(C), ADcholC)

x0 <- seq(0,1,length=1000)

Bx <- LMMsolver:::Bsplines(knots, x0)
v1 <- sigma2e*(diag(Bx %*% solve(C) %*% t(Bx)))
v2 <- diag(Bx %*% sparseInverse %*% t(Bx))

v1
v2
range(v1-v2)

X <- cbind(1,x)
U <- cbind(X,B)
UtU <- crossprod(U)
constr1 <- c(1,rep(0, q-1))
constr2 <- c(rep(0,q-1),1)
constr <- cbind(constr1,constr2)
C <- (UtU + lambda*bdiag.spam(diag(0,2), DtD + tcrossprod(constr) ))/sigma2e

X0 <- cbind(1, x0)
U0 <- cbind(X0, Bx)

v3 <- diag(U0 %*% solve(C) %*% t(U0))

lC <- list()
lC[[1]] <- UtU
lC[[2]] <- bdiag.spam(diag(0,2), DtD + tcrossprod(constr) )
ADcholC <- LMMsolver:::ADchol(lC)
ED <- theta * LMMsolver:::dlogdet(ADcholC, theta)

sparseInverse2 <- LMMsolver:::DerivCholesky(chol(C), ADcholC)

v4 <- diag(U0 %*% sparseInverse2 %*% t(U0))

v5 <- rowSums((U0 %*% obj$sparseInverse) * U0)

range(v1-v2)
range(v1-v3)
range(v1-v4)
range(v1-v5)

plotDat$se <- sqrt(v5)
plotDat$ysim <- simFun(x0)

ggplot(data = dat, aes(x = x, y = y)) +
  geom_point(size = 1.2) +
  geom_line(data = plotDat, aes(y = ypred), color = "red", size = 1) +
  geom_line(data = plotDat, aes(y = ysim), color = "green", size = 1) +
  geom_line(data = plotDat, aes(y = ypred-2*se), col='blue', size=1) +
  geom_line(data = plotDat, aes(y = ypred+2*se), col='blue', size=1) +
  theme(panel.grid = element_blank())
