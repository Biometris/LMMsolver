load("testdata.rda")

## Fit models with different components.
mod <- LMMsolve(fixed = y ~ 1,
                spline = ~spl3D(x1 = x1, x2 = x2, x3 = x3,
                                nseg = c(4, 4, 4)),
                data = simDat,
                tolerance = 1e-3)

## Check summary function.
dims <- summary(mod, which = "dimensions")
expect_inherits(dims, "data.frame")
expect_equivalent_to_reference(dims, "effDims")

varcomps <- summary(mod, which = "variances")
expect_inherits(varcomps, "data.frame")
expect_equivalent_to_reference(varcomps, "varComps")

expect_stdout(print(summary(mod)),
              "Table with effective dimensions and penalties")
expect_stdout(print(summary(mod)),
              "Total Effective Dimension: 250")

expect_stdout(print(summary(mod, which = "variances")),
              "Table with variances")

## Check logLik function.
expect_equal(logLik(mod), 198.972874670064)
expect_equal(logLik(mod, includeConstant = FALSE), 421.355999705594)

## Check deviance function.
expect_equal(deviance(mod), -397.945749340128)
expect_equal(deviance(mod, includeConstant = FALSE), -842.71199941116)

## Check coef function.
coefMod <- coef(mod)
expect_inherits(coefMod, "list")
expect_equal(names(coefMod), c("(Intercept)", "splF", "splR"))

expect_equal_to_reference(coefMod, "modCoefs", tolerance = 1e-6)

## Check fitted function.
fitMod <- fitted(mod)

expect_equal_to_reference(fitMod, "modFit", tolerance = 1e-6)

## Check residuals function.
residMod <- residuals(mod)

expect_equal_to_reference(residMod, "modResid", tolerance = 1e-6)