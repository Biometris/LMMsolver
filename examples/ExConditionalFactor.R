library(asreml)
library(dplyr)
library(LMMsolver)
data(john.alpha, package = "agridat")
head(john.alpha)
dat <- john.alpha
Nblk <- nlevels(dat$block)*nlevels(dat$rep)
dat$block2 <- as.factor(rep(paste0("BLK",  formatC(1:Nblk,width=2,flag="0")),each=4))

obj0 <- asreml(fixed = yield ~ gen+rep,
               random = ~block:at(rep),
               data = dat, tolerance=1.0e-10)
summary(obj0)

obj1 <- asreml(fixed = yield ~ gen+rep,
               random = ~at(rep,"R1"):block + at(rep,"R2"):block,
               data = dat, tolerance=1.0e-10)
summary(obj1)

# by foot in LMMsolver
dat_ext <- dat
dat_ext$R1_cond <- 1*(dat$rep=="R1")
dat_ext$R2_cond <- 1*(dat$rep=="R2")

obj2 <- LMMsolve(fixed=yield~gen+rep,
                random=~block:R1_cond + block:R2_cond,data=dat_ext)
summary(obj2)
obj1$loglik
obj2$logL

obj3 <- LMMsolve(fixed=yield~gen+rep,
                 random=~cf(block,rep,"R1") + cf(block,rep,"R2"),
                 data=dat)
summary(obj3)

obj1$loglik
obj2$logL
obj3$logL

