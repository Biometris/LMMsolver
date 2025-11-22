## Use the oats data described in JABES2920 paper for some tests of spl1D.
data("oats.data")
load("testdata.rda")

## Fit models as described in JABES2020 paper.

## Baseline model, only fixed effects.
obj0 <- LMMsolve(fixed = yield ~ rep + gen,
                 data = oats.data,
                 tolerance = 1.0e-10)

## Number of plots
N <- nrow(oats.data)

## Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield ~ rep + gen,
                 spline = ~spl1D(x = plot, nseg = N - 1, degree = 1, pord = 1,
                                 scaleX = FALSE),
                 data = oats.data,
                 tolerance = 1.0e-10)

## Test input checks in obtainSmoothTrend.
expect_error(obtainSmoothTrend("obj0"),
             "object should be an object of class LMMsolve")
expect_error(obtainSmoothTrend(obj0),
             "The model was fitted without a spline component")
expect_error(obtainSmoothTrend(obj1),
             "Specify either grid or newdata")
expect_error(obtainSmoothTrend(obj1, grid = 72, deriv = -1),
             "deriv should be an integer greater than or equal to zero")
expect_error(obtainSmoothTrend(obj1, grid = 72, deriv = 2),
             "2-order derivatives cannot be computed for B-splines of")

## Trend using grid.

expect_error(obtainSmoothTrend(obj1, grid = c(1, 2)),
             "grid should be a numeric vector with length equal to the")

obj1Trend1 <- obtainSmoothTrend(obj1, grid = 72)

## Trend using newdata.

expect_error(obtainSmoothTrend(obj1, newdata = "oats.data"),
             "newdata should be a data.frame")
expect_error(obtainSmoothTrend(obj1, newdata = oats.data[, -1]),
             "The following smoothing variables are not in newdata")

obj1Trend2 <- obtainSmoothTrend(obj1, newdata = oats.data)

## Results should be the same.
expect_equal(obj1Trend1$ypred, obj1Trend2$ypred)

## Include intercept.
obj1Trend3 <- obtainSmoothTrend(obj1, newdata = oats.data,
                                includeIntercept = TRUE)
expect_equivalent(obj1Trend3$ypred - obj1Trend2$ypred,
                  rep(coef(obj1)$`(Intercept)`, N))

obj2 <- LMMsolve(fixed = yield ~ 1, random = ~ gen,
                 spline = ~ spl2D(x1 = bed, x2 = row, nseg = c(3, 3)),
                 data = durban.rowcol,
                 tolerance = 1e-6)

## Trend using grid.

expect_error(obtainSmoothTrend(obj2, grid = 2),
             "grid should be a numeric vector with length equal to the")
expect_warning(obtainSmoothTrend(obj2, grid = c(5, 5), deriv = 1),
               "deriv is ignored for 2-dimensional splines")

expect_error(obtainSmoothTrend(obj2, grid = c(5, 5), which = 3),
        "which should be an integer with value at most the number of fittedspline components.\n")

obj2Trend1 <- obtainSmoothTrend(obj2, grid = c(5, 5))

## Trend using newdata.

obj2Trend2 <- obtainSmoothTrend(obj2, newdata = durban.rowcol)
## Reorder output to match order of predictions on grid.
obj2Trend2 <- obj2Trend2[order(obj2Trend2$bed, obj2Trend2$row), ]

## Results should be the same.
expect_equal(obj2Trend1$ypred, obj2Trend2$ypred)


## Use simulated data for testing spl3D.
load("testdata.rda")

obj3 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3, nseg = c(4, 4, 4)),
                 data = simDat,
                 tolerance = 1e-3)

## Trend using grid.

expect_error(obtainSmoothTrend(obj3, grid = 2),
             "grid should be a numeric vector with length equal to the")
expect_warning(obtainSmoothTrend(obj3, grid = c(5, 5, 5), deriv = 1),
               "deriv is ignored for 3-dimensional splines")

obj3Trend1 <- obtainSmoothTrend(obj3, grid = c(5, 5, 5))

expect_equal_to_reference(obj3Trend1, "smooth3D1")

## Trend using newdata.

obj3Trend2 <- obtainSmoothTrend(obj3, newdata = simDat)

expect_equal_to_reference(obj3Trend2, "smooth3D2")


