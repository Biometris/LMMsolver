# predict example...
rm(list=ls())

library(agridat)
library(LMMsolver)
library(ggplot2)
library(asreml)
library(spam)

data(john.alpha)
dat <- john.alpha
# remove one observation, to make example unbalanced:
dat <- dat[-20,]

obj1 <- asreml(fixed = yield~rep+gen,
                 random=~rep:block,
                 data = dat)

obj2 <- LMMsolve(fixed = yield~rep+gen,
                random=~rep:block,
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
summary(obj2)

# asreml
coef(obj1, list=TRUE)$'rep:block'
# labels missing
coef(obj2)$'rep:block'

displayMME(obj2)

obj1$loglik
obj2$logL

# predictions using asreml
pred1 <- predict(obj1, classify='gen',sed=TRUE)
#head(pred1$pvals[1:4,])
pred1$sed[1:3,1:3]

# predictions plus standard errors,
# averaging over fixed rep, ignoring random rep:block
# to make this general, we have to store fixed/random info in LMMsolve object
Nrep <- nlevels(dat$rep)
Ngen <- nlevels(dat$gen)
Nblk <- nlevels(dat$block)*Nrep
Dg1 <- matrix(c(1, rep(1/Nrep, Nrep-1), 0, 0, rep(0, Ngen-3), rep(0, Nblk)), nrow=1)
Dg2 <- matrix(c(1, rep(1/Nrep, Nrep-1), 1, 0, rep(0, Ngen-3), rep(0, Nblk)), nrow=1)
Dg3 <- matrix(c(1, rep(1/Nrep, Nrep-1), 0, 1, rep(0, Ngen-3), rep(0, Nblk)), nrow=1)

Dg <- rbind.spam(Dg1, Dg2, Dg3)
pred2se <- LMMsolver:::calcStandardErrors(obj2$C, Dg)
pred1se <- pred1$pvals$std.error[1:3]
range(pred1se - pred2se)

# sed
Dg12 <- Dg1-Dg2
Dg13 <- Dg1-Dg3
Dg23 <- Dg2-Dg3
Dg2 <- rbind.spam(Dg12, Dg13, Dg23)
pred2sed <- LMMsolver:::calcStandardErrors(obj2$C, Dg2)

A <- matrix(data=NA,ncol=3,nrow=3)
A[1,2] <- A[2,1] <- pred2sed[1]
A[1,3] <- A[3,1] <- pred2sed[2]
A[2,3] <- A[3,2] <- pred2sed[3]
A

range(A - pred1$sed[1:3,1:3], na.rm=TRUE)

