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

dat1 <- dat[1:36, ]
dat2 <- dat[37:72,]

# splines by foot:
txtSpl <- as.character(spl)[2]
splRes1 <- eval(parse(text = txtSpl), envir = dat1, enclos = parent.frame())
q1 <- ncol(splRes1$Z)

# splines by foot:
splRes2 <- eval(parse(text = txtSpl), envir = dat2, enclos = parent.frame())
q2 <- ncol(splRes2$Z)

lGrp <- list(B1 = 1:q, B2=c((q+1):(2*q)))
lGinv <- list(B1 = splRes1$lGinv[[1]], B2=splRes2$lGinv[[1]])

x1 <- splRes1$X[ , 1]
x2 <- splRes2$X[ , 1]

zero1 <- matrix(data=0,ncol=1,nrow=36)
zero2 <- matrix(data=0,ncol=q1,nrow=36)

p1 <- rbind(as.matrix(splRes1$Z), zero2)
p2 <- rbind(zero2, as.matrix(splRes2$Z))
x1_ext <- c(x1, zero1)
x2_ext <- c(zero1, x2)

dat_ext <- cbind(p1,p2,x1_ext,x2_ext, dat)

obj2 <- LMMsolve(fixed=yield~ gen+rep + x1_ext + x2_ext,
                 random = ~grp(B1)+grp(B2),
                 ginverse = lGinv,
                 group = lGrp,
                 data=dat_ext)

summary(obj2)

dat$half <- as.factor(rep(c("H1","H2"),each=36))

obj3 <- LMMsolve(fixed=yield~gen+rep,
                 spline=~spl1D(row,nseg=nseg,cond=half, level="H1") +
                         spl1D(row,nseg=nseg,cond=half, level="H2"),
                 data=dat)

summary(obj3)

obj2$logL
obj3$logL

