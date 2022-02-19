# Simulate data, using first example JOPS book:
library(ggplot2)
library(LMMsolver)
library(spam)
library(agridat)

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
obj <-LMMsolve(y~1, spline=~spl1D(x,nseg=nseg, scaleX=FALSE, xlim=c(0,1)),data=dat)
summary(obj)

# some calculations by foot, using the parameters theta obtained from LMMsolve
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
partialDerivChol <- LMMsolver:::DerivCholesky(chol(C), ADcholC)

x0 <- seq(0,1,length=1000)
Bx <- LMMsolver:::Bsplines(knots, x0)
v1 <- sigma2e*(diag(Bx %*% solve(C) %*% t(Bx)))
v2 <- diag(Bx %*% partialDerivChol %*% t(Bx))

# calculate mixed model matrix C for mixed model sparse formulation:
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
partialDerivChol2 <- LMMsolver:::DerivCholesky(chol(C), ADcholC)

v4 <- diag(U0 %*% partialDerivChol2 %*% t(U0))
v5 <- rowSums((U0 %*% obj$partialDerivChol) * U0)

v7 <- sapply(1:nrow(U0), FUN= function(x) {sum(partialDerivChol2 * crossprod(U0[x,]))})
range(v4-v7)

plotDat <- obtainSmoothTrend(obj, grid = 1000, includeIntercept = TRUE)
head(plotDat)
v6 <- (plotDat$se)^2

eig <- eigen(DtD)
Usc <- eig$vectors[, -c(q-1,q)] %*% diag(eig$values[-c(q-1,q)]^(-1/2))
Z <- B %*% Usc
dat_ext <- cbind(dat, Z)

L <- list(ndx = c(3:18))
obj2 <- LMMsolve(y~x,random=~grp(ndx), group=L, data=dat_ext)

obj$logL - obj2$logL

U <- cbind(X, Z)
UtU <- crossprod(U)
C <- (UtU + lambda*bdiag.spam(diag(0,2), diag(1,q-2) ))/sigma2e

X0 <- cbind(1, x0)
U0 <- cbind(X0, Bx %*% Usc)

v7 <- diag(U0 %*% solve(C) %*% t(U0))

range(v1-v2)
range(v1-v3)
range(v1-v4)
range(v1-v5)
range(v1-v6)
range(v1-v7)

plotDat$ysim <- simFun(x0)

ggplot(data = dat, aes(x = x, y = y)) +
  geom_point(size = 1.2) +
  geom_line(data = plotDat, aes(y = ypred), color = "red", size = 1) +
  geom_line(data = plotDat, aes(y = ysim), color = "green", size = 1) +
  geom_line(data = plotDat, aes(y = ypred-2*se), col='blue', size=1) +
  geom_line(data = plotDat, aes(y = ypred+2*se), col='blue', size=1) +
  theme(panel.grid = element_blank())


# other examples, with more fixed (or random terms) terms
data(john.alpha)
dat <- john.alpha

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield~rep,
                 random = ~gen,
                 spline = ~spl1D(x = plot, nseg = nseg, scaleX=FALSE),
                 data = dat,
                 trace = FALSE,
                 tolerance = 1.0e-10)
summary(obj1)

# obtain spatial trend with genotype fixed:
plotDat <- obtainSmoothTrend(obj1, grid=100)
head(plotDat)

head(plotDat)
ggplot(plotDat, aes(x = plot, y = ypred)) +
  geom_line() +
  geom_line(data = plotDat, aes(y = ypred-2*se), col='blue', size=1) +
  geom_line(data = plotDat, aes(y = ypred+2*se), col='blue', size=1) +
  labs(title = "Spatial trend for the oats data", x = "plot", y = "plot effects") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


obj1$term.labels.f
obj1$term.labels.r
labels <- c(obj1$term.labels.f, obj1$term.labels.r)
obj1$dim

splRes <- obj1$splRes[[1]]

e <- cumsum(obj1$dim)
s <- e - obj1$dim + 1

# if splRes$term != NULL, otherwise skip...
ndx.f <- which(splRes$term.labels.f == labels)
ndx.r <- which(splRes$term.labels.r == labels)

# choose grid....
x0 <- seq(1,72,by=3)
Bx <- LMMsolver:::Bsplines(splRes$knots[[1]], x0)

U <- spam(x=0, ncol=sum(obj1$dim),nrow=length(x0))
U[,1] <- 1
U[, c(s[ndx.f]:e[ndx.f])] <- x0
U[, c(s[ndx.r]:e[ndx.r])] <- Bx
display(U)

