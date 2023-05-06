# test for using same penalty for variance components:

library(LMMsolver)
library(spam)
set.seed(1234)

# number of observations:
n <- 100
# number of levels for three components:
nL <- c(25, 20, 15)

# sd for rnorm
sd <- c(2, 2, 1)

a_eff <- rnorm(nL[1], sd=sd[1])
b_eff <- rnorm(nL[2], sd=sd[2])
c_eff <- rnorm(nL[3], sd=sd[3])

sA <- sample(1:nL[1], size=n, replace=TRUE)
sB <- sample(1:nL[2], size=n, replace=TRUE)
sC <- sample(1:nL[3], size=n, replace=TRUE)

eps <- rnorm(n, 1.0)
y <- a_eff[sA] + b_eff[sB] + c_eff[sC] + eps

dat <- data.frame(y=y, A=as.factor(sA), B=as.factor(sB), C=as.factor(sC))
head(dat)

obj0 <- LMMsolve(y~1, random=~A+B+C, data=dat)
summary(obj0)

obj1 <- LMMsolve(y~1, random=~A+B+C, data=dat, grpTheta=c(1,4,2,3))
summary(obj1)

# should give same result
obj0$logL - obj1$logL

# give A and B the same penalty:
obj2 <- LMMsolve(y~1, random=~A+B+C, data=dat, grpTheta=c(1,1,2,3))
summary(obj2)

# give A, B, C all the same penalty:
obj3 <- LMMsolve(y~1, random=~A+B+C, data=dat, grpTheta=c(1,1,1,2))
summary(obj3)



