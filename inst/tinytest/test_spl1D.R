load("testdata.rda")

# use the oats data described in JABES2920 paper for some tests of spl1D.

# baseline model, only fixed effects.
obj0 <- LMMsolve(fixed = yield ~ rep + gen,
                 data = john.alpha,
                 trace = FALSE,
                 tolerance = 1.0e-10,
                 omitConstant = FALSE)

## compare with JABES2020 paper, table 1:
devJABES2020paper_baseline <- 69.91
expect_equal(round(obj0$dev, 2), devJABES2020paper_baseline)

# number of plots
N <- nrow(john.alpha)

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield ~ rep + gen,
                 spline = ~spl1D(x = plot, nseg = N - 1, degree = 1, pord = 1,
                                 scaleX = FALSE),
                 data = john.alpha,
                 trace = FALSE,
                 tolerance = 1.0e-10,
                 omitConstant = FALSE)

## compare with JABES2020 paper LV model, table 1:
devJABES2020paper_LV <- 54.49
expect_equal(round(obj1$dev,2) , devJABES2020paper_LV)
