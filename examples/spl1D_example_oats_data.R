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

obj1 <- LMMsolve(fixed = yield~rep+gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
obj1$ED

# calculate deviance, as in JABES 2020 paper, table 1:
p <- 1 + (v-1) + (r-1)
Constant = log(2*pi)*(N-p)
dev = -2*obj$logL + Constant
round(dev,2)

# residual variance, see JABES 2020 paper, table 1:
round(obj$sigma2e, 5)

# genotype random, not in JABES paper
obj2 <- LMMsolve(fixed = yield~rep,
                 random = ~gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
obj2$ED
