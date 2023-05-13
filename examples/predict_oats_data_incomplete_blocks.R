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

whatPred = "gen"

# predictions using asreml
pred1 <- predict(obj1, classify=whatPred, sed=TRUE)
#head(pred1$pvals[1:4,])
pred1$sed[1:3,1:3]

pred2 <- LMMsolver:::predictTest(obj2, classify=whatPred)

range(pred1$pvals$predicted.value - pred2$prediction)
range(pred1$pvals$std.error - pred2$se)

obj3 <- asreml(fixed = yield~rep,
               random=~gen+rep:block,
               data = dat)

obj4 <- LMMsolve(fixed = yield~rep,
                 random=~gen+rep:block,
                 data = dat,
                 trace = FALSE,
                 tolerance = 1.0e-10)
summary(obj4)

obj3$loglik
obj4$logL

# predictions using asreml
pred3 <- predict(obj3, classify=whatPred, sed=TRUE)
#head(pred1$pvals[1:4,])
pred3$sed[1:3,1:3]

pred4 <- LMMsolver:::predictTest(obj4, classify=whatPred)

range(pred3$pvals$predicted.value - pred4$prediction)
range(pred3$pvals$std.error - pred4$se)

# predictions plus standard errors,
# averaging over fixed rep, ignoring random rep:block
# to make this general, we have to store fixed/random info in LMMsolve object
Nrep <- nlevels(dat$rep)
Ngen <- nlevels(dat$gen)
Nblk <- nlevels(dat$block)*Nrep

# matrix for geno fixed, column 1 is reference:
Mgen <- diag(Ngen)[,-1]
Mmu  <- matrix(data=1, nrow=Ngen, ncol=1)
# averaging over fixed effect replicate:
Mrep <- matrix(data=rep(1/Nrep, Nrep-1), nrow=Ngen, ncol=Nrep-1, byrow=TRUE)
# ignoring random effect rep:block
Mblk <- matrix(data=0, nrow=Ngen, ncol=Nblk)
# make the matrix for predictions genotype
Dg <- cbind.spam(Mmu, Mrep, Mgen, Mblk)

pred1se <- pred1$pvals$std.error
pred2se <- LMMsolver:::calcStandardErrors(obj2$C, Dg)
range(pred1se - pred2se)

# make the predictors for sed's:
DgSED_all <- NULL
for (i in 2:Ngen) {
  for (j in 1:i) {
    if (i!=j) {
      DgSED <- Dg[i, ] - Dg[j, ]
      DgSED_all <- rbind.spam(DgSED_all, DgSED)
    }
  }
}

# just to see the structure of VCOV:
tmp <- obj2$C + 0 * crossprod(DgSED_all)
display(tmp)

pred2sed <- LMMsolver:::calcStandardErrors(obj2$C, DgSED_all)
A <- matrix(data=NA, nrow=Ngen, ncol=Ngen)
k = 1
for (i in 2:Ngen) {
  for (j in 1:i) {
    if (i!=j) {
      A[i,j] <- pred2sed[k]
      k <- k+1
    }
  }
}
range(A - pred1$sed, na.rm = TRUE)

