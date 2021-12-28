load("testdata.rda")

## Some tests on correct input checking for spl3D.

## Test that formula is specified in a correct way.
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = yield ~ spl3D(x1 = x1, x2 = x2, x3 = x3,
                                             nseg = c(4, 4, 4))),
             "spline should be a formula of form")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3,
                                      nseg = c(4, 4, 4)) +
                        spl4D(x1 = x1, x2 = x2, x3 = x3, nseg = c(20, 20, 20))),
             "spline should be a formula of form")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = tst, x2 = x2, x3 = x3,
                                      nseg = c(4, 4, 4))),
             "The following variables in the spline part of the model are not")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = tst,
                                      nseg = c(4, 4, 4))),
             "The following variables in the spline part of the model are not")

## Test other arguments of spl3D.
# nseg
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3,
                                      nseg = c(-1, 10, 10))),
             "nseg should be a vector of length three of positive integers")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = 10)),
             "nseg should be a vector of length three of positive integers")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3,
                                      nseg = c(2.5, 10, 10))),
             "nseg should be a vector of length three of positive integers")

# pord
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      pord = 3)),
             "pord should be either 1 or 2")

# degree
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      degree = -1)),
             "degree should be a positive integer")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      degree = c(2, 2))),
             "degree should be a positive integer")

# x1lim
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x1lim = 10)),
             "x1lim should be a vector of length two")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x1lim = c(2, 72))),
             "x1lim should be a vector of length two")

# x2lim
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x2lim = 10)),
             "x2lim should be a vector of length two")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x2lim = c(2, 72))),
             "x2lim should be a vector of length two")

# x3lim
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x3lim = 10)),
             "x3lim should be a vector of length two")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x3lim = c(2, 72))),
             "x3lim should be a vector of length two")

## Fit model.

obj1 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4)),
                 data = simDat,
                 tolerance = 1e-3)

## Check that full LMM solve object is correct.
expect_equivalent_to_reference(obj1, "spl3DFull")
