## Use the oats data described in JABES2920 paper for some tests of spl1D.
data("john.alpha", package = "agridat")

## Some tests on correct input checking for spl1D.

## Test that formula is specified in a correct way.
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = yield ~ spl1D(x = plot, nseg = 10)),
             "spline should be a formula of form")
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = 10) +
                        spl4D(x = plot, nseg = 20)),
             "spline should be a formula of form")
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = tst, nseg = 10) ),
             "The following variables in the spline part of the model are not")

## Test other arguments of spl1D.
# nseg
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = -1)),
             "nseg should be a positive integer")
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = c(2, 2))),
             "nseg should be a positive integer")

# pord
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = 10, pord = 3)),
             "pord should be either 1 or 2")

# degree
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = 10, degree = -1)),
             "degree should be a positive integer")
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = 10, degree = c(2, 2))),
             "degree should be a positive integer")

# xlim
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = 10, xlim = 10)),
             "xlim should be a vector of length two")
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = 10, xlim = c(2, 72))),
             "xlim should be a vector of length two")

## Fit models as described in JABES2020 paper.

## Baseline model, only fixed effects.
obj0 <- LMMsolve(fixed = yield ~ rep + gen,
                 data = john.alpha,
                 tolerance = 1.0e-10)

## Compare deviance with JABES2020 paper, table 1.
devJABES2020paper_baseline <- 69.91
expect_equal(round(deviance(obj0), 2), devJABES2020paper_baseline)

## Number of plots
N <- nrow(john.alpha)

## Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield ~ rep + gen,
                 spline = ~spl1D(x = plot, nseg = N - 1, degree = 1, pord = 1,
                                 scaleX = FALSE),
                 data = john.alpha,
                 tolerance = 1.0e-10)

## Compare deviance with JABES2020 paper LV model, table 1.
devJABES2020paper_LV <- 54.49
expect_equal(round(deviance(obj1), 2) , devJABES2020paper_LV)

## Check that full LMM solve object is correct.
expect_equivalent_to_reference(obj1, "spl1DFull")

