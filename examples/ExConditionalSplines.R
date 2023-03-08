library(dplyr)
library(LMMsolver)
library(spam)
data(john.alpha, package = "agridat")
head(john.alpha)
dat <- john.alpha

nseg <- 50
degr <- 3

spl <- ~spl1D(row,nseg=nseg)
obj0 <- LMMsolve(fixed=yield~gen+rep, spline=spl, data=dat)

# splines by foot:
txtSpl <- as.character(spl)[2]
splRes <- eval(parse(text = txtSpl), envir = dat, enclos = parent.frame())
q <- ncol(splRes$Z)
lGrp <- list(B = 1:q)
lGinv <- list(B = splRes$lGinv[[1]])
x <- splRes$X[ , 1]

# combine the elements in extended data frame
dat_ext <- cbind(as.matrix(splRes$Z), x, dat)

obj1 <- LMMsolve(fixed=yield~ gen + rep + x,
                 random = ~grp(B),
                 ginverse = lGinv,
                 group = lGrp,
                 data=dat_ext)

# results and comparisons:
summary(obj0)
summary(obj1)

obj0$logL
obj1$logL

#
# splines within replicate
#

# by foot, using group and ginverse arguments
dat1 <- dat[dat$rep=='R1', ]
dat2 <- dat[dat$rep=='R2', ]
dat3 <- dat[dat$rep=='R3', ]

txtSpl <- as.character(spl)[2]
splRes1 <- eval(parse(text = txtSpl), envir = dat1, enclos = parent.frame())
q1 <- ncol(splRes1$Z)

splRes2 <- eval(parse(text = txtSpl), envir = dat2, enclos = parent.frame())
q2 <- ncol(splRes2$Z)

splRes3 <- eval(parse(text = txtSpl), envir = dat3, enclos = parent.frame())
q3 <- ncol(splRes3$Z)

lGrp <- list(B1 = 1:q, B2=c((q+1):(2*q)), B3=c((2*q+1):(3*q)))
lGinv <- list(B1 = splRes1$lGinv[[1]],
              B2 = splRes2$lGinv[[1]],
              B3 = splRes3$lGinv[[1]])

x1 <- splRes1$X[ , 1]
x2 <- splRes2$X[ , 1]
x3 <- splRes3$X[ , 1]

zeroV <- rep(1, 24)
zeroM <- matrix(data=0,ncol=q1,nrow=24)

p1 <- rbind(as.matrix(splRes1$Z), zeroM, zeroM)
p2 <- rbind(zeroM, as.matrix(splRes2$Z),zeroM)
p3 <- rbind(zeroM, zeroM, as.matrix(splRes3$Z))

x1_ext <- c(x1, zeroV, zeroV)
x2_ext <- c(zeroV, x2, zeroV)
x3_ext <- c(zeroV, zeroV, x3)

dat_ext <- cbind(p1,p2,p3,x1_ext,x2_ext,x3_ext, dat)

obj2 <- LMMsolve(fixed=yield~ gen + rep + x1_ext + x2_ext + x3_ext,
                 random = ~grp(B1)+grp(B2) + grp(B3),
                 ginverse = lGinv,
                 group = lGrp,
                 data=dat_ext)

summary(obj2)

# splines within replicate using conditional formatting splines:
obj3 <- LMMsolve(fixed=yield~gen+rep,
                 spline=~spl1D(row,nseg=nseg, cond=rep, level="R1") +
                         spl1D(row,nseg=nseg, cond=rep, level="R2") +
                         spl1D(row,nseg=nseg, cond=rep, level="R3"),
                 data=dat)

summary(obj3)

obj2$logL
obj3$logL

