## Use spring oats data with lattice design.
data(oats.data)

## Checks for correct inputs.
obj <- LMMsolve(yield ~ rep,
                      random = ~gen + block:gen,
                      data = oats.data)

expect_error(displayMME(object = obj$C),
             "object should be an object of class LMMsolve.")

expect_error(diagnosticsMME(object = obj$C),
             "object should be an object of class LMMsolve.")
