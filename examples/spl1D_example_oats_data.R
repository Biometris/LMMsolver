# Martin Boer, Biometris, Wageningen
# example sparse 1D, see Boer et al 2020, JABES for details
# of the example
rm(list=ls())

library(agridat)
library(LMMsolver)
library(ggplot2)

data(john.alpha)
dat <- john.alpha

N <- nrow(dat)

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield~rep+gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1, scaleX=FALSE),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
summary(obj1)

# use LV model....
cN <- c(1/sqrt(N-1), rep(0, N-2),1/sqrt(N-1))
D <- diff(diag(N), diff = 1)
Delta <- 0.5*crossprod(D)
LVinv <- 2*Delta + cN %*% t(cN)
lGinv <- list(plotF=LVinv)

dat$plotF <- as.factor(dat$plot)

obj1b <- LMMsolve(fixed = yield~rep+gen,
                 random = ~plotF,
                 ginverse = lGinv,
                 data = dat,
                 trace = FALSE,
                 tolerance = 1.0e-10)
summary(obj1b)

dev1  <- deviance(obj1)
dev1b <- deviance(obj1b)
dev1 - dev1b

# residual variance, see JABES 2020 paper, table 1:
round(obj1$sigma2e, 5)

# genotype random, not in JABES paper
obj2 <- LMMsolve(fixed = yield~rep,
                 random = ~gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1, scaleX=FALSE),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
summary(obj2)

# genotype random, not in JABES paper
obj2b <- LMMsolve(fixed = yield~rep,
                 random = ~gen+plotF,
                 ginverse = lGinv,
                 data = dat,
                 trace = FALSE,
                 tolerance = 1.0e-10)
# genotype random, not in JABES paper
obj2c <- LMMsolve(fixed = yield~rep,
                  random = ~plotF+gen,
                  ginverse = lGinv,
                  data = dat,
                  trace = FALSE,
                  tolerance = 1.0e-10)

summary(obj2b)

dev2  <- deviance(obj2)
dev2b <- deviance(obj2b)
dev2c <- deviance(obj2c)
dev2 - dev2b
dev2 - dev2c

diagnosticsMME(obj2)
displayMME(obj2, cholesky=FALSE)
displayMME(obj2, cholesky=TRUE)

displayMME(obj2b, cholesky=FALSE)


# obtain spatial trend with genotype fixed:
plotDat <- obtainSmoothTrend(obj1, grid=100)
head(plotDat)
ggplot(plotDat, aes(x = plot, y = ypred)) +
  geom_line() +
  labs(title = "Spatial trend for the oats data", x = "plot", y = "plot effects") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


## Obtain special trend by specifying newdat
newdat <- data.frame(plot = seq(1, 72, length.out = 100))

plotDat2 <- obtainSmoothTrend(obj1, newdata = newdat)

## Results should be identical
identical(plotDat$ypred, plotDat2$ypred)

ggplot(plotDat2, aes(x = plot, y = ypred)) +
  geom_line() +
  labs(title = "Spatial trend for the oats data", x = "plot", y = "plot effects") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

