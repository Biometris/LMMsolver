#
# modelling splines within replicates.
#
library(dplyr)
library(LMMsolver)
library(spam)
data(john.alpha, package = "agridat")
head(john.alpha)
dat <- john.alpha

# number of segments per replicate, here for testing
# chosen to be different
nseg <- c(20, 15, 10)
degr <- 3
q <- nseg + 3

# by foot, using group and ginverse arguments
dat1 <- dat[dat$rep=='R1', ]
dat2 <- dat[dat$rep=='R2', ]
dat3 <- dat[dat$rep=='R3', ]

txtSpl1 <- paste0("spl1D(row, nseg = ", nseg[1], ")")
txtSpl2 <- paste0("spl1D(row, nseg = ", nseg[2], ")")
txtSpl3 <- paste0("spl1D(row, nseg = ", nseg[3], ")")

splRes1 <- eval(parse(text = txtSpl1), envir = dat1)
splRes2 <- eval(parse(text = txtSpl2), envir = dat2)
splRes3 <- eval(parse(text = txtSpl3), envir = dat3)

lGrp <- list(B1 = c(1:q[1]),
             B2 = c((q[1]+1):(q[1]+q[2])),
             B3 = c((q[1]+q[2]+1):(q[1]+q[2]+q[3])))

lGinv <- list(B1 = splRes1$lGinv[[1]],
              B2 = splRes2$lGinv[[1]],
              B3 = splRes3$lGinv[[1]])

x1 <- splRes1$X[ , 1]
x2 <- splRes2$X[ , 1]
x3 <- splRes3$X[ , 1]
zeroV <- rep(1, 24)
x1_ext <- c(x1, zeroV, zeroV)
x2_ext <- c(zeroV, x2, zeroV)
x3_ext <- c(zeroV, zeroV, x3)

Z <- matrix(data=0, ncol=sum(q), nrow=72)
Z[1:24,  c(1:q[1])] <- as.matrix(splRes1$Z)
Z[25:48, c((q[1]+1):(q[1]+q[2]))] <- as.matrix(splRes2$Z)
Z[49:72, c((q[1]+q[2]+1):(sum(q)))] <- as.matrix(splRes3$Z)

dat_ext <- cbind(Z, x1_ext,x2_ext,x3_ext, dat)

obj2 <- LMMsolve(fixed=yield~ gen + rep + x1_ext + x2_ext + x3_ext,
                 random = ~grp(B1)+grp(B2) + grp(B3),
                 ginverse = lGinv,
                 group = lGrp,
                 data=dat_ext)

summary(obj2)

# splines within replicate using conditional formatting splines:
obj3 <- LMMsolve(fixed=yield~gen+rep,
                 spline=~spl1D(row,nseg=nseg[1], cond=rep, level="R1") +
                         spl1D(row,nseg=nseg[2], cond=rep, level="R2") +
                         spl1D(row,nseg=nseg[3], cond=rep, level="R3") ,
                 data=dat)

summary(obj3)

obj2$logL
obj3$logL

obj2$logL - obj3$logL

