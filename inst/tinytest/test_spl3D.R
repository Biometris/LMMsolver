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
                                      pord = 4)),
             "pord should be equal to 1, 2 or 3")

expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      degree = 2, pord = 3)),
             "pord should be less or equal to degree.\n")

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
             "All values of x1 should be between the lower and upper value of x1lim")

# x2lim
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x2lim = 10)),
             "x2lim should be a vector of length two")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x2lim = c(2, 72))),
             "All values of x2 should be between the lower and upper value of x2lim")

# x3lim
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x3lim = 10)),
             "x3lim should be a vector of length two")
expect_error(LMMsolve(fixed = y ~ 1, data = simDat,
                      spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4),
                                      x3lim = c(2, 72))),
             "All values of x3 should be between the lower and upper value of x3lim")

## Fit model.

obj1 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4)),
                 data = simDat,
                 tolerance = 1e-3)

## From R 4.3 there is an extra item in the family output.
## This gives problems with the comparison.
## Therefore it is removed first.

obj1$family$dispersion <- NULL

## Check that full LMM solve object is correct.
expect_equivalent_to_reference(obj1, "spl3DFull")

# pord = 1, no fixed effect
obj2 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4), pord = 1),
                 data = simDat,
                 tolerance = 1e-3)
expect_equivalent_to_reference(obj2, "spl3DFull2")

# pord = 3
obj3 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4), pord = 3),
                 data = simDat,
                 tolerance = 1e-3)
expect_equivalent_to_reference(obj3, "spl3DFull3", tolerance = 1e-6)
