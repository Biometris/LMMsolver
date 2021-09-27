library(agridat)

data(john.alpha)
dat <- john.alpha

n <- 4     # number of units per block
b <- 6     # number of blocks per replicate
r <- 3     # number of replicates
v <- 24    # number of genotypes/replicate:
N <- n*b*r # total number of observations.

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper:
obj1 <- LMMsolve(fixed = yield~rep+gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1, scaleX=FALSE),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
obj1$ED

# calculate deviance, as in JABES 2020 paper, table 1,
# include the constant omitted in LMMsolver (and asreml, SpATS)
p <- 1 + (v-1) + (r-1)
Constant = log(2*pi)*(N-p)
dev = -2*obj1$logL + Constant

# compare with JABES2020 paper
devJABES2020paper <- 54.49
expect_equal(round(dev,2), devJABES2020paper)

