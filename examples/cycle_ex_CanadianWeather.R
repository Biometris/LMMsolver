#' ---
#' title: Circular data Montreal Temperature
#' author: Martin Boer, Biometris, WUR, Wageningen
#' ---

library(ggplot2)
library(LMMsolver)
suppressMessages(library(spam))
suppressMessages(library(fda))

data(CanadianWeather)

df <- as.data.frame(CanadianWeather$dailyAv[,,1])
dim(df)

city <- "Montreal"

x <- c(1:365)/365
y <- df[[city]]

dat <- data.frame(day=x, temp=y)

nseg <- 100
tol <- 1.0e-6

# Prepare basis and penalty
knots <- LMMsolver:::PsplinesKnots(0, 1, degree = 3, nseg = nseg, cyclic = TRUE)
cB <- LMMsolver:::Bsplines(knots, x)
q <- ncol(cB)
cD <- LMMsolver:::cDiff(q)
cD[1:3,1:3]

G <- LMMsolver:::constructG(knots,scaleX = FALSE, pord=2)
range(cD %*% G[,-1])
X <- cB %*% G
x1 <- X[,2]
x2 <- X[,3]
range(x1^2 + x2^2-1)
plot(x=x1, y=x2)

#' Method 1: Only fixed part (null-space) of model
#' ==================
#'
dat1 <- data.frame(y, x1, x2)
obj1 <- LMMsolve(fixed = y~x1+x2, data = dat1)
summary(obj1)

# make predictions on a grid
x0 <- seq(0, 1, by=0.001)
cB0 <- LMMsolver:::Bsplines(knots,x0)
X0 <- cB0 %*% G
X0 <- X0[,-1]

# make predictions by foot:
a1 <- obj1$coefMME
Z1 <- cbind(1, X0)
eta1 <- as.vector(Z1 %*% a1)
mu1 <- eta1

#' Method 2: Eigen decomposition of penalty
#' ==================

# calculate eigen decomposition of penalty P
P <- t(cD) %*% cD
eig <- eigen(P)
d_plus <- eig$values[-c(q-1,q)]
U_plus <- eig$vectors[,-c(q-1,q)]
U_sc <- U_plus %*% diag(1/sqrt(d_plus))
BU <- cB %*% U_sc
dat2 <- data.frame(BU, y, x1, x2)

# Including the random (non-linear part), by defining the
# columns in the dataframe dat
lGrp <- list(ndx = c(1:ncol(BU)))
obj2 <- LMMsolve(fixed = y~x1+x2,
                 random = ~grp(ndx), # using grp()
                 group = lGrp,       # defines the list
                 data = dat2,
                 tolerance = tol)
summary(obj2)

a2 <- obj2$coefMME
BU0 <- cB0 %*% U_sc
Z2 <- cbind(1, X0, BU0)
eta2 <- as.vector(Z2 %*% a2)
mu2 <- eta2

#' Method 3: Add constraints to random effect
#' ==================

dat3 <- data.frame(as.matrix(cB), y, x1, x2)
q <- ncol(cB)
# add additional constraints, see Boer (2023).
cB_star <- LMMsolver:::Bsplines(knots, x=c(0,0.25))

#P <- as.spam(crossprod(cD) + crossprod(cB_star))
scaleFactor <- LMMsolver:::calcScaleFactor(list(knots), pord=2)
P <- LMMsolver:::constructGinvSplines(q,list(knots),pord=2,scaleFactor=scaleFactor)
lGinv <- list(ndx = P[[1]])
lGrp <- list(ndx = c(1:q))
obj3 <- LMMsolve(fixed = y~x1 +x2,
                 random = ~grp(ndx), # using grp()
                 group = lGrp,       # defines the list
                 ginverse = lGinv,   # the penalty matrix
                 data = dat3,
                 tolerance = tol,
                 trace=TRUE)
summary(obj3)

# the sum to zero built-in constraints for random effect
cB_star %*% coef(obj3)$ndx

a3 <- obj3$coefMME
Z3 <- cbind(1, X0, cB0)
eta3 <- as.vector(Z3 %*% a3)
mu3 <- eta3

# should give same results
range(mu2-mu3)


#' Method 4: Use spl1D with cyclic option
#' ==================

dat4 <- data.frame(y, x)

obj4 <- LMMsolve(fixed = y~1,
                 spline = ~spl1D(x=x,nseg=nseg, cyclic=TRUE, scaleX=FALSE),
                 data = dat4,
                 tolerance = tol, trace=TRUE)
summary(obj4)

# compare models 2 and 3 and 4:
deviance(obj2)
deviance(obj3)
deviance(obj4)

newdat <- data.frame(x=x0)
pred <- predict(obj4, newdata=newdat)
range(mu3-pred$ypred)

#' Make plots
#' ==================

Data <- data.frame(x = c(x0,x0), mu=c(mu1,mu3),
                   model=as.factor(c(rep("fixed", length(x0)),
                         rep("circ P-splines",length(x0)))))

plt2 <- ggplot(Data, aes(x = x, y=mu, col=model))+
  ggtitle(paste(city, ": January thaw (Ramsay, Silverman 2005, p73)")) +
  xlab("year") + ylab("temperature") +
  geom_line() +
  geom_point(data=dat,aes(x=day,y=temp), inherit.aes = FALSE,size=0.5)

print(plt2)

# some extra code....
spl <- obj4$splRes[[1]]
all.equal(spl$X, X[,-1])
all.equal(spl$Z, cB)
obj2$nIter
obj3$nIter
obj4$nIter



