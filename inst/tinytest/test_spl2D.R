## Use the oats data described in BioRxiv 2021 paper by Hans-Peter Piepho.
data("durban.rowcol", package = "agridat")
## Restict data to 5 columns and rows for speed of testing.
durban.rowcol <- durban.rowcol[durban.rowcol$row < 6 & durban.rowcol$bed < 6, ]

## Some tests on correct input checking for spl2D.

## Test that formula is specified in a correct way.
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = yield ~ spl2D(x1 = bed, x2 = row, nseg = c(10, 10))),
             "spline should be a formula of form")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10)) +
                        spl4D(x1 = bed, x2 = row, nseg = c(20, 20))),
             "spline should be a formula of form")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = tst, x2 = row, nseg = c(10, 10))),
             "The following variables in the spline part of the model are not")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = tst, nseg = c(10, 10))),
             "The following variables in the spline part of the model are not")

## Test other arguments of spl2D.
# nseg
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(-1, 10))),
             "nseg should be a vector of length two of positive integers")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = 10)),
             "nseg should be a vector of length two of positive integers")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(2.5, 10))),
             "nseg should be a vector of length two of positive integers")

# pord
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      pord = 3)),
             "pord should be either 1 or 2")

# degree
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      degree = -1)),
             "degree should be a positive integer")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      degree = c(2, 2))),
             "degree should be a positive integer")

# x1lim
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x1lim = 10)),
             "x1lim should be a vector of length two")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x1lim = c(2, 72))),
             "x1lim should be a vector of length two")

# x2lim
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x2lim = 10)),
             "x2lim should be a vector of length two")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x2lim = c(2, 72))),
             "x2lim should be a vector of length two")

## Fit simplified (for speed) version of model described in bioRxiv 2021 paper.

obj1 <- LMMsolve(fixed = yield ~ 1, random = ~ gen,
                 spline = ~ spl2D(x1 = bed, x2 = row, nseg = c(3, 3)),
                 data = durban.rowcol,
                 tolerance = 1e-6)

## Check that full LMM solve object is correct.
expect_equivalent_to_reference(obj1, "spl2DFull", tolerance = 1e-6)
