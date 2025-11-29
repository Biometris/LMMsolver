## Fit model on oats data
data(oats.data)

obj <- LMMsolve(fixed  = yield ~ rep,
                random = ~gen,
                spline = ~spl1D(x = plot, nseg = 20),
                data = oats.data)
EDdf <- effDim(obj)

expect_equivalent_to_reference(EDdf, "effDim0")

# example from Schmidt et al., Genetics 2019:
obj2 <- LMMsolve(fixed = yield ~ rep,
                  random = ~gen + rep:block,
                  data = oats.data)
EDdf2 <- effDim(obj2)
expect_equal(round(subset(EDdf2, Term == "gen")$Ratio, 3), 0.809)

