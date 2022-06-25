library(LMMsolver)
library(asreml)
library(tictoc)
library(spam)

set.seed(1234)

# total number of observations
n <- 250000

# number of levels for the two variance components
nL <- c(50, 100, 150)

G1 <-as.factor(sample(1:nL[1], size=n, replace=TRUE))
G2 <-as.factor(sample(1:nL[2], size=n, replace=TRUE))
G3 <-as.factor(sample(1:nL[3], size=n, replace=TRUE))

dat <- data.frame(G1=G1, G2=G2, G3=G3)
Z1 <- model.matrix(~G1-1,dat)
Z2 <- model.matrix(~G2-1,dat)
Z3 <- model.matrix(~G3-1,dat)

g1_eff <- rnorm(n=nL[1], sd= 1.0)
g2_eff <- rnorm(n=nL[2], sd=10.0)
g3_eff <- rnorm(n=nL[3], sd=25.0)

e <- rnorm(n, sd=100.0)
y <- Z1 %*% g1_eff + Z2 %*% g2_eff + Z3 %*% g3_eff + e
dat$y <- y

tic()
obj0 <- LMMsolve(fixed=y~1,random=~G1+G2+G3, data=dat, trace=TRUE)
toc()

summary(obj0)
summary(obj0,which='variances')

tic()
obj1 <- asreml(fixed=y~1, random=~G1+G2+G3, data=dat, trace=TRUE)
toc()

obj0$logL
obj1$loglik

obj0$logL - obj1$loglik
