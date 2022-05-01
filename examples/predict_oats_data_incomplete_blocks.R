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

# Not  efficient way to calculate SEDs, just to compare results
# with asreml
A <- matrix(data=NA, nrow=Ngen, ncol=Ngen)
for (i in 1:Ngen) {
  for (j in 1:Ngen) {
    if (i!=j) {
      DgSED <- Dg[i,] - Dg[j, ]
      pred2sed <- LMMsolver:::calcStandardErrors(obj2$C, DgSED)
      A[i,j] <- pred2sed[1]
    }
  }
}
range(A - pred1$sed, na.rm = TRUE)

