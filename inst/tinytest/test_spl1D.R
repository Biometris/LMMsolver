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
                      spline = ~spl1D(x = plot, nseg = 10, xlim = "a")),
             "xlim should be numeric")
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = 10, xlim = 10)),
             "xlim should be a vector of length two")
expect_error(LMMsolve(fixed = yield ~ rep + gen, data = john.alpha,
                      spline = ~spl1D(x = plot, nseg = 10, xlim = c(2, 72))),
             "All values of plot should be between the lower and upper value of xlim")

## Fit models as described in JABES2020 paper.

## Baseline model, only fixed effects.
obj0 <- LMMsolve(fixed = yield ~ rep + gen,
                 data = john.alpha,
                 tolerance = 1.0e-10)

## Compare deviance with JABES2020 paper, table 1.
devJABES2020paper_baseline <- 69.91
expect_equal(round(deviance(obj0, relative=FALSE), 2), devJABES2020paper_baseline)

## Number of plots
N <- nrow(john.alpha)

## Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield ~ rep + gen,
                 spline = ~spl1D(x = plot, nseg = N - 1, degree = 1, pord = 1,
                                 scaleX = FALSE),
                 data = john.alpha, tolerance = 1e-5)

## Compare deviance with JABES2020 paper LV model, table 1.
devJABES2020paper_LV <- 54.49
expect_equal(round(deviance(obj1, relative=FALSE), 2) , devJABES2020paper_LV)

## From R 4.3 there is an extra item in the family output.
## This gives problems with the comparison.
## Therefore it is removed first.

obj1$family$dispersion <- NULL

## Check that full LMM solve object is correct.
expect_equivalent_to_reference(obj1, "spl1DFull")

## Test that LMMsolver:: can be used inside spline part.
expect_silent(LMMsolve(fixed = yield ~ rep + gen,
                       spline = ~LMMsolver::spl1D(x = plot, nseg = N - 1),
                       data = john.alpha))

## Test combination of splines with conditional factor.
obj2 <- LMMsolve(fixed = yield ~ rep + gen,
                 spline= ~spl1D(row, nseg = 25, cond = rep, level = "R1") +
                   spl1D(row, nseg = 20, cond = rep, level = "R2") +
                   spl1D(row, nseg = 15, cond = rep, level = "R3"),
                 data = john.alpha)

sumObj2 <- summary(obj2)
expect_equal(nrow(sumObj2), 10)
expect_equal(sumObj2[["Term"]],
             c("(Intercept)", "rep", "gen", "lin(row)_R1", "lin(row)_R2",
               "lin(row)_R3", "s(row)_R1", "s(row)_R2", "s(row)_R3", "residual"))
expect_equal(sumObj2[["Ratio"]],
             c(1, 1, 1, 1, 1, 1, 0.325761468207249, 0.0865062712483912,
               0.0593138565407634, 0.738710428505654))



