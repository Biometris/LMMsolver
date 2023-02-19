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

obj1$loglik
obj2$logL

# compare coefficient, se and zRatio for genotype
asr <- summary(obj1, coef = TRUE)
cf <- coef(obj2, se=TRUE)
asr.gen <- asr$coef.fixed[1:24,]

range(cf$gen$value  - asr.gen[,1])
range(cf$gen$se     - asr.gen[,2], na.rm = TRUE)
range(cf$gen$zRatio - asr.gen[,3], na.rm = TRUE)
