load("testdata.rda")

## Test that general input checks work correctly.
expect_error(LMMsolve(fixed = pheno ~ cross, data = "testDat"),
             "data should be a data.frame")

## Test fixed, random and residual formulas (spline part is tested in spl checks).
expect_error(LMMsolve(fixed = ~ cross, data = testDat),
             "fixed should be a formula of the form")
expect_error(LMMsolve(fixed = pheno ~ cross, random = pheno ~ cross,
                      data = testDat),
             "random should be a formula of the form")
expect_error(LMMsolve(fixed = pheno ~ cross, residual = pheno ~ cross,
                      data = testDat),
             "residual should be a formula of the form")

expect_error(LMMsolve(fixed = pheno ~ tst, data = testDat),
             "The following variables in the fixed part of the model are not")
expect_error(LMMsolve(fixed = pheno ~ cross, random = ~tst, data = testDat),
             "The following variables in the random part of the model are not")
expect_error(LMMsolve(fixed = pheno ~ cross, residual = ~tst, data = testDat),
             "The following variables in the residual part of the model are not")

## Test ginverse.
ginv <- matrix(1:4, nrow = 2)
ginvL <- list(ginv = ginv)
ginvLS <- list(ginv = ginv %*% t(ginv))
expect_error(LMMsolve(fixed = pheno ~ cross, ginverse = ginv, data = testDat),
             "ginverse should be a named list of symmetric matrices")
expect_error(LMMsolve(fixed = pheno ~ cross, ginverse = ginvL, data = testDat),
             "ginverse should be a named list of symmetric matrices")
expect_error(LMMsolve(fixed = pheno ~ cross, ginverse = ginvLS, data = testDat),
             "ginverse element ginv not defined in random part")

## Test other input parameters.
expect_error(LMMsolve(fixed = pheno ~ cross, data = testDat, tolerance = -1),
             "tolerance should be a positive numerical value")
expect_error(LMMsolve(fixed = pheno ~ cross, data = testDat, maxit= -1),
             "maxit should be a positive numerical value")

## Test use of grp - group.
Lgrp <- list(QTL = 3:5)
expect_error(LMMsolve(fixed = pheno ~ cross, group = Lgrp, data = testDat),
             "The following variables in group are not specified in grp")
expect_error(LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL2),
                      data = testDat),
             "The following variables are specified in grp in the random part")
expect_error(LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL2),
                      group = Lgrp, data = testDat),
             "The following variables are specified in grp in the random part")
expect_error(LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL),
                      group = c(Lgrp, list(QTL2 = 1:2)), data = testDat),
             "The following variables in group are not specified in grp")

## Fit models with different components.
mod0 <- LMMsolve(fixed = pheno ~ cross, data = testDat)
mod1 <- LMMsolve(fixed = pheno ~ cross, residual = ~cross, data = testDat)
mod2 <- LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL),
                 group = Lgrp, data = testDat)
mod3 <- LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL) + ind,
                 group = Lgrp, data = testDat)
mod4 <- LMMsolve(fixed = pheno ~ cross, random = ~ind, data = testDat)

## Compare results with predefined output.
expect_equivalent_to_reference(mod0, "LMMsolve0")
expect_equivalent_to_reference(mod1, "LMMsolve1")
expect_equivalent_to_reference(mod2, "LMMsolve2")
expect_equivalent_to_reference(mod3, "LMMsolve3")
expect_equivalent_to_reference(mod4, "LMMsolve4")

