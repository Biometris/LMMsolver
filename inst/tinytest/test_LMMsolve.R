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
ginvLS2 <- list(ind = ginv %*% t(ginv))
indMat <- diag(nrow = nlevels(testDat$ind))
rownames(indMat) <- colnames(indMat) <- levels(testDat$ind)
ginvLS3 <- list(ind = indMat)
expect_error(LMMsolve(fixed = pheno ~ cross, ginverse = ginv, data = testDat),
             "ginverse should be a named list of symmetric matrices")
expect_error(LMMsolve(fixed = pheno ~ cross, ginverse = ginvL, data = testDat),
             "ginverse should be a named list of symmetric matrices")
expect_error(LMMsolve(fixed = pheno ~ cross, ginverse = ginvLS, data = testDat),
             "ginverse element ginv not defined in random part")
expect_error(LMMsolve(fixed = pheno ~ cross, random = ~ind, ginverse = ginvLS2,
                      data = testDat),
             "Dimensions of ind should match number of levels")

## Test other input parameters.
expect_error(LMMsolve(fixed = pheno ~ cross, data = testDat, tolerance = -1),
             "tolerance should be a positive numerical value")
expect_error(LMMsolve(fixed = pheno ~ cross, data = testDat, maxit= -1),
             "maxit should be a positive numerical value")
expect_warning(LMMsolve(fixed = pheno ~ cross, data = testDat, maxit= 1),
               "No convergence after 1 iterations")

## Test use of grp - group.
Lgrp <- list(QTL = 3:5)
expect_error(LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL2),
                      data = testDat),
             "The following variables are specified in grp in the random part")
expect_error(LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL2),
                      group = Lgrp, data = testDat),
             "The following variables are specified in grp in the random part")

## Fit models with different components.
mod0 <- LMMsolve(fixed = pheno ~ cross, data = testDat)
mod1 <- LMMsolve(fixed = pheno ~ cross, residual = ~cross, data = testDat)
mod2 <- LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL),
                 group = Lgrp, data = testDat)
mod3 <- LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL) + ind,
                 group = Lgrp, data = testDat)
mod4 <- LMMsolve(fixed = pheno ~ cross, random = ~ind, data = testDat)
mod5 <- LMMsolve(fixed = pheno ~ cross, random = ~ind, ginverse = ginvLS3,
                 data = testDat)

## Compare results with predefined output.
expect_equivalent_to_reference(mod0, "LMMsolve0")
expect_equivalent_to_reference(mod1, "LMMsolve1")
expect_equivalent_to_reference(mod2, "LMMsolve2")
expect_equivalent_to_reference(mod3, "LMMsolve3")
expect_equivalent_to_reference(mod4, "LMMsolve4")
expect_equivalent_to_reference(mod5, "LMMsolve5")

## Group not used in random part should be ignored.
mod3a <- LMMsolve(fixed = pheno ~ cross, random = ~grp(QTL) + ind,
                  group = c(Lgrp, list(QTL2 = 1:2)), data = testDat)
expect_equivalent(mod3, mod3a)

## Test option trace.
expect_stdout(LMMsolve(fixed = pheno ~ cross, data = testDat, trace = TRUE),
              "iter logLik")
expect_stdout(LMMsolve(fixed = pheno ~ cross, data = testDat, trace = TRUE),
              "2 -207.1218")

## Test that zero-variance in response in caught.
testDatZv1 <- testDatZv2 <- testDat
testDatZv1[["pheno"]] <- 1
testDatZv2[testDatZv2[["cross"]] == "AxB", "pheno"] <- 1

expect_error(LMMsolve(fixed = pheno ~ cross, data = testDatZv1),
             "Variance response variable zero or almost zero")
expect_error(LMMsolve(fixed = pheno ~ cross, residual = ~cross,
                      data = testDatZv2),
             "Variance response variable zero or almost zero for levels")

## Test that variables with only NA are caught.
testDatNA <- testDat
testDatNA[["cross"]] <- NA
expect_error(LMMsolve(fixed = pheno ~ cross, data = testDatNA),
             "in the fixed part of the model only have missing values")

## Test that result for character variables is identical to that of factors.
testDatChar <- testDat
testDatChar[["cross"]] <- as.character(testDatChar[["cross"]])

modChar <-  LMMsolve(fixed = pheno ~ cross, residual = ~cross,
                     data = testDatChar)

## Checking equality gives a warning since spam doesn't check attributes.
expect_equivalent(modChar, mod1)

## Test that interaction terms are labeled correctly.
testDatInt <- testDat
testDatInt[["rep"]] <- factor(c(rep(c("r1", "r2"), each = 50),
                                rep(c("r1", "r2"), each = 40)))

modInt <- LMMsolve(fixed = pheno ~ rep:cross, data = testDatInt)
coefsInt <- coef(modInt)$`rep:cross`
expect_equal(names(coefsInt),
             c("rep_r1:cross_AxB", "rep_r2:cross_AxB", "rep_r1:cross_AxC",
               "rep_r2:cross_AxC"))
expect_equivalent(coefsInt,
                  c(-0.438185444999992, -0.158384204999993, 0.394857450000009, 0))

expect_silent(LMMsolve(fixed = pheno ~ rep + rep:cross, data = testDatInt))
expect_silent(LMMsolve(fixed = pheno ~ rep*cross, data = testDatInt))

## Test that interaction terms are labeled correctly when level missing.
testDatInt2 <- testDat
testDatInt2[["rep"]] <- factor(c(rep(c("r1", "r2"), each = 50),
                                 rep("r1", times = 80)))

# Interaction in fixed part.
modIntMiss <- LMMsolve(fixed = pheno ~ rep:cross, data = testDatInt2)
coefsIntMiss <- coef(modIntMiss)$`rep:cross`
expect_equal(names(coefsIntMiss),
             c("rep_r1:cross_AxB", "rep_r2:cross_AxB", "rep_r1:cross_AxC"))
expect_equivalent(coefsIntMiss,
                  c(-0.63561417, -0.355812930000001, 0))

# Interaction in random part.

modIntMissR <- LMMsolve(fixed = pheno ~ 1, random = ~rep:cross,
                        data = testDatInt2)
coefsIntRMiss <- coef(modIntMissR)$`rep:cross`
expect_equal(names(coefsIntRMiss),
             c("rep_r1:cross_AxB", "rep_r2:cross_AxB", "rep_r1:cross_AxC",
               "rep_r2:cross_AxC"))
expect_equivalent(coefsIntRMiss,
                  c(-0.136895043027739, -0.0221165394103219,
                    0.159011582438064, 0))


