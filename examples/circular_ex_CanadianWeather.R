#' ---
#' title: Circular data Montreal Temperature
#' author: Martin Boer, Biometris, WUR, Wageningen
#' ---

library(ggplot2)
library(LMMsolver)
suppressMessages(library(spam))
suppressMessages(library(fda))

data(CanadianWeather)

df <- as.data.frame(CanadianWeather$dailyAv[,,1])
dim(df)

city <- "Montreal"

x <- c(1:365)/365
y <- df[[city]]

dat <- data.frame(day=x, temp=y)

nseg <- 100

# Prepare basis and penalty
knots <- LMMsolver:::PsplinesKnots(0, 1, degree = 3, nseg = nseg)
cB <- LMMsolver:::cBsplines(knots, x)
q <- ncol(cB)
cD <- LMMsolver:::cDiff(q)

# define null - space
z <- c(0:(nseg-1))/nseg
u1 <- sin(2*pi*z)
u2 <- cos(2*pi*z)
range(cD %*% u1)
range(cD %*% u2)
x1 <- cB %*% u1
x2 <- cB %*% u2

#' Method 1: Only fixed part (null-space) of model
#' ==================
#'
dat1 <- data.frame(y, x1, x2)
obj1 <- LMMsolve(fixed = y~x1+x2, data = dat1)
summary(obj1)

# make predictions on a grid
x0 <- seq(0, 1, by=0.001)
cB0 <- LMMsolver:::cBsplines(knots,x0)
x1_grid <- cB0 %*% u1
x2_grid <- cB0 %*% u2

# make predictions by foot:
a1 <- obj1$coefMME
Z1 <- cbind(1, x1_grid, x2_grid)
eta1 <- as.vector(Z1 %*% a1)
mu1 <- eta1

#' Method 2: Eigen decomposition of penalty
#' ==================

# calculate eigen decomposition of penalty P
P <- t(cD) %*% cD
eig <- eigen(P)
d_plus <- eig$values[-c(q-1,q)]
U_plus <- eig$vectors[,-c(q-1,q)]
U_sc <- U_plus %*% diag(1/sqrt(d_plus))
BU <- cB %*% U_sc
dat2 <- data.frame(BU, y, x1, x2)

# Including the random (non-linear part), by defining the
# columns in the dataframe dat
lGrp <- list(ndx = c(1:ncol(BU)))
obj2 <- LMMsolve(fixed = y~x1+x2,
                 random = ~grp(ndx), # using grp()
                 group = lGrp,       # defines the list
                 data = dat2)
summary(obj2)

a2 <- obj2$coefMME
BU0 <- cB0 %*% U_sc
Z2 <- cbind(1, x1_grid, x2_grid, BU0)
eta2 <- as.vector(Z2 %*% a2)
mu2 <- eta2

#' Method 3: Add constraints to random effect
#' ==================

dat3 <- data.frame(as.matrix(cB), y, x1, x2)
q <- ncol(cB)
# add additional constraints, see Boer (2023).
cB_star <- LMMsolver:::cBsplines(knots, x=c(0,0.25))

P <- as.spam(crossprod(cD) + crossprod(cB_star))
lGinv <- list(ndx = P)
lGrp <- list(ndx = c(1:q))
obj3 <- LMMsolve(fixed = y~x1 +x2,
                 random = ~grp(ndx), # using grp()
                 group = lGrp,       # defines the list
                 ginverse = lGinv,   # the penalty matrix
                 data = dat3)
summary(obj3)

# the sum to zero built-in constraints for random effect
cB_star %*% coef(obj3)$ndx

# compare models 2 and 3:
deviance(obj2)
deviance(obj3)

a3 <- obj3$coefMME
Z3 <- cbind(1, x1_grid, x2_grid, cB0)
eta3 <- as.vector(Z3 %*% a3)
mu3 <- eta3

# should give same results
range(mu2-mu3)

#' Make plots
#' ==================

Data <- data.frame(x = c(x0,x0), mu=c(mu1,mu3),
                   model=as.factor(c(rep("fixed", length(x0)),
                         rep("circ P-splines",length(x0)))))

plt2 <- ggplot(Data, aes(x = x, y=mu, col=model))+
  ggtitle(paste(city, ": January thaw (Ramsay, Silverman 2005, p73)")) +
  xlab("year") + ylab("temperature") +
  geom_line() +
  geom_point(data=dat,aes(x=day,y=temp), inherit.aes = FALSE,size=0.5)

print(plt2)

