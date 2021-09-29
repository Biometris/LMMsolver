# Martin Boer, Biometris, Wageningen
# example sparse 1D, see Boer et al 2020, JABES for details
# of the example
rm(list=ls())

library(agridat)
library(LMMsolver)

data(john.alpha)
dat <- john.alpha

n <- 4     # number of units per block
b <- 6     # number of blocks per replicate
r <- 3     # number of replicates
v <- 24    # number of genotypes/replicate:
N <- n*b*r # total number of observations.

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield~rep+gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1, scaleX=FALSE),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10,
                omitConstant = FALSE)
obj1$ED

round(obj1$dev, 2)

# residual variance, see JABES 2020 paper, table 1:
round(obj1$sigma2e, 5)

# genotype random, not in JABES paper
obj2 <- LMMsolve(fixed = yield~rep,
                 random = ~gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1, scaleX=FALSE),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
obj2$ED

tr <- obtainSmoothTrend1D(obj2,grid=100)
str(tr)
plot(x=tr[[1]]$x, y=tr$eta,type='l')


