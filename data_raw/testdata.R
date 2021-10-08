## Simulate data for testing 3D splines.
set.seed(1234)

Wood3D <- function(x1, x2, x3) { # Wood 2006
  y <- 1.5 * exp(-(x1 - 0.2) ^ 2 / 5 - (x2-0.5) ^ 2 / 3 - (x3 - 0.9) ^ 2 / 4) +
    0.5 * exp(-(x1 - 0.3) ^ 2 / 4 -(x2 - 0.7) ^ 2 / 2 - (x3 - 0.4) ^ 2 / 6) +
    exp(-(x1 - 0.1) ^ 2 / 5 -(x2 - 0.3) ^ 2 / 5 - (x3 - 0.7) ^ 2 / 4)
  return(y)
}

N <- 250
x1 <- runif(N)
x2 <- runif(N)
x3 <- runif(N)
eps <- rnorm(N, sd = 0.1)

y <- Wood3D(x1, x2, x3) + eps

simDat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

## Save to tinytest directory.
save(simDat, file = "./inst/tinytest/testdata.rda")
