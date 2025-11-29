## Fit model on oats data
data(oats.data)

obj <- LMMsolve(fixed  = yield ~ rep,
                random = ~gen,
                spline = ~spl1D(x = plot, nseg = 20),
                data = oats.data)
EDdf <- effDim(obj)

expect_equivalent_to_reference(EDdf, "effDim0")

