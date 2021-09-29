load("testdata.rda")

# use the oats data described in JABES2920 paper for some tests of spl1D.

n <- 4         # number of units per block
b <- 6         # number of blocks per replicate
r <- 3         # number of replicates
v <- 24        # number of genotypes/replicate:
N <- n * b * r # total number of observations.

## include the constant for logL omitted in LMMsolver (and asreml, SpATS).
p <- 1 + (v - 1) + (r - 1)
Constant = log(2 * pi) * (N - p)

# baseline model, only fixed effects.
obj0 <- LMMsolve(fixed = yield ~ rep + gen,
                 data = john.alpha,
                 trace = FALSE,
                 tolerance = 1.0e-10)
dev0 = -2 * obj0$logL + Constant

## compare with JABES2020 paper, table 1:
devJABES2020paper_baseline <- 69.91
expect_equal(round(dev0, 2), devJABES2020paper_baseline)

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield ~ rep + gen,
                 spline = ~spl1D(x = plot, nseg = N - 1, degree = 1, pord = 1,
                                 scaleX = FALSE),
                 data = john.alpha,
                 trace = FALSE,
                 tolerance = 1.0e-10)
dev1 = -2 * obj1$logL + Constant

## compare with JABES2020 paper LV model, table 1:
devJABES2020paper_LV <- 54.49
expect_equal(round(dev1, 2), devJABES2020paper_LV)
