## Use spring oats data with lattice design.
data(oats.data)

## Checks for correct inputs.
expect_error(LMMsolve(yield ~ gen,
                      random = ~ cf(var = blocky, cond = rep, level = "R1"),
                      data = oats.data),
             "blocky should be a variable in data")

expect_error(LMMsolve(yield ~ gen,
                      random = ~ cf(var = block, cond = repy, level = "R1"),
                      data = oats.data),
             "repy should be a variable in data")

expect_error(LMMsolve(yield ~ gen,
                      random = ~ cf(var = block, cond = plot, level = "R1"),
                      data = oats.data),
             "cond should be a factor")

expect_error(LMMsolve(yield ~ gen,
                      random = ~ cf(var = block, cond = rep, level = NULL),
                      data = oats.data),
             "level should be defined")

expect_error(LMMsolve(yield ~ gen,
                      random = ~ cf(var = block, cond = rep, level = "R4"),
                      data = oats.data),
             "R4 is not a level of rep")


obj0 <- LMMsolve(fixed = yield ~ gen + rep,
                 random = ~ cf(var = block, cond = rep, level = "R1") +
                   cf(var = block, cond = rep, level = "R2"),
                 data = oats.data)

## Check that full LMM solve object is correct.

## From R 4.3 there is an extra item in the family output.
## This gives problems with the comparison.
## Therefore it is removed first.

obj0$family$dispersion <- NULL

expect_equivalent_to_reference(obj0, "cfFull")

## Check that NA in response doesn't crash the model (issue #102)
ja <- oats.data
ja[1, "yield"] <- NA

expect_warning(LMMsolve(fixed = yield ~ gen + rep,
                        random = ~ cf(var = block, cond = rep, level = "R1") +
                          cf(var = block, cond = rep, level = "R2"),
                        data = ja),
               "1 observations removed with missing value for yield")





