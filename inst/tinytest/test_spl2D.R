
# part of durban.rowcol experiment, first 5 columns and rows for speed of testing.
durban.rowcol <- structure(list(row = c(4L, 5L, 3L, 4L, 4L, 4L, 1L, 3L, 2L, 2L,
                                        3L, 5L, 3L, 1L, 3L, 4L, 5L, 2L, 1L, 2L, 1L, 1L, 5L, 2L, 5L),
                                bed = c(4L, 2L, 4L, 5L, 3L, 1L, 5L, 2L, 1L, 4L, 3L, 5L, 5L,
                                        4L, 1L, 2L, 1L, 5L, 3L, 2L, 1L, 2L, 4L, 3L, 3L),
                                rep = structure(c(1L,
                                                  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), levels = c("R1", "R2"
                                                  ), class = "factor"),
                                gen = structure(1:25, levels = c("G002", "G005", "G007",
                                                                 "G008", "G012", "G019", "G039", "G042", "G050", "G099", "G107",
                                                                 "G125", "G133", "G137", "G144", "G168", "G190", "G191", "G192",
                                                                 "G194", "G198", "G201", "G232", "G252", "G270"), class = "factor"),
                                yield = c(5.28, 5.14, 5.77, 4.94, 5.49, 5.54, 5.1, 5.86,
                                          5.46, 4.83, 6.1, 4.72, 5.88, 5.23, 5.51, 5.87, 5.38, 5.8,
                                          6.07, 5.81, 5.83, 5.78, 4.88, 5.12, 5.35)),
                           row.names = c(1:25), class = "data.frame")

## Some tests on correct input checking for spl2D.

## Test that formula is specified in a correct way.
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = yield ~ spl2D(x1 = bed, x2 = row, nseg = c(10, 10))),
             "spline should be a formula of form")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10)) +
                        spl4D(x1 = bed, x2 = row, nseg = c(20, 20))),
             "spline should be a formula of form")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = tst, x2 = row, nseg = c(10, 10))),
             "The following variables in the spline part of the model are not")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = tst, nseg = c(10, 10))),
             "The following variables in the spline part of the model are not")

## Test other arguments of spl2D.
# nseg
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(-1, 10))),
             "nseg should be a vector of length two of positive integers")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = 10)),
             "nseg should be a vector of length two of positive integers")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(2.5, 10))),
             "nseg should be a vector of length two of positive integers")

# pord
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      pord = 4)),
             "pord should be equal to 1, 2 or 3")

# pord
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      pord = 3, degree = 2)),
             "pord should be less or equal to degree")


# degree
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      degree = -1)),
             "degree should be a positive integer")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      degree = c(2, 2))),
             "degree should be a positive integer")

# x1lim
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x1lim = 10)),
             "x1lim should be a vector of length two")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x1lim = c(2, 72))),
             "All values of bed should be between the lower and upper value of x1lim")

# x2lim
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x2lim = 10)),
             "x2lim should be a vector of length two")
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x2lim = c(2, 72))),
             "All values of row should be between the lower and upper value of x2lim")

# cyclic
expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x2lim = c(2, 72), cyclic = c(3,4))),
             "cyclic should be a logical vector of length two", fixed=TRUE)

expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x2lim = c(2, 72), cyclic= c(TRUE,FALSE))),
             "x1 should be in the range [0,1] for cyclic data", fixed=TRUE)

expect_error(LMMsolve(fixed = yield ~ 1, data = durban.rowcol,
                      spline = ~spl2D(x1 = bed, x2 = row, nseg = c(10, 10),
                                      x2lim = c(2, 72), cyclic= c(FALSE,TRUE))),
             "x2 should be in the range [0,1] for cyclic data", fixed=TRUE)


## Fit simplified (for speed) version of model described in bioRxiv 2021 paper.

obj1 <- LMMsolve(fixed = yield ~ 1, random = ~ gen,
                 spline = ~ spl2D(x1 = bed, x2 = row, nseg = c(3, 3)),
                 data = durban.rowcol,
                 tolerance = 1e-6)

## From R 4.3 there is an extra item in the family output.
## This gives problems with the comparison.
## Therefore it is removed first.

obj1$family$dispersion <- NULL

## Check that full LMM solve object is correct.
expect_equivalent_to_reference(obj1, "spl2DFull", tolerance = 1e-6)

## Test combination of splines with conditional factor.

obj2 <- LMMsolve(fixed = yield ~ 1, random = ~ gen,
                 spline = ~ spl2D(x1 = bed, x2 = row, nseg = c(3, 3),
                                  cond = rep, level = "R1"),
                 data = durban.rowcol,
                 tolerance = 1e-6)

sumObj2 <- summary(obj2)
expect_equal(nrow(sumObj2), 6)
expect_equal(sumObj2[["Term"]],
             c("(Intercept)", "lin(bed, row)_R1", "gen", "s(bed)_R1",
               "s(row)_R1", "residual"))
expect_equal(sumObj2[["Ratio"]],
             c(1, 1, 0.458895940871542, 4.82648592298555e-05,
               0.0821598533833015, 0.458895940885926))

## pord=1
obj1a <- LMMsolve(fixed = yield ~ 1, random = ~ gen,
                 spline = ~ spl2D(x1 = bed, x2 = row, nseg = c(3, 3), pord = 1),
                 data = durban.rowcol,
                 tolerance = 1e-6)
expect_equivalent_to_reference(obj1a, "spl2DFull2", tolerance= 1.0e-6)

## pord=3
obj1b <- LMMsolve(fixed = yield ~ 1, random = ~ gen,
                  spline = ~ spl2D(x1 = bed, x2 = row, nseg = c(3, 3) , pord = 3),
                  data = durban.rowcol,
                  tolerance = 1e-6)
expect_equivalent_to_reference(obj1b, "spl2DFull3", tolerance= 1.0e-6)

## cyclic: cylinder

set.seed(1234)

sim_fun <- function(x) {
  y <- exp(sin(1.25*pi*x[1])) + exp(sin(2*pi*(x[2]-0.1)))
  y
}

# simulate data on a grid:
n1 <- 100
n2 <- 100
n <- n1*n2
dat <- expand.grid(x1=seq(0,2,length=n1),x2=seq(0,1,length=n2))
dat$ytrue <- apply(dat, MARGIN=1, FUN = sim_fun)

# take subset for training
Ntraining <- 2500
k <- sample(x=c(1:n), size = Ntraining)
dat_train <- dat[k, ]
dat_train$y <- dat_train$ytrue + rnorm(n=Ntraining)

nseg <- c(20,20)

obj3 <- LMMsolve(fixed = y~1,
                spline = ~spl2D(x1=x1,x2=x2, nseg=nseg,
                                cyclic=c(FALSE,TRUE)),
                data = dat_train)
expect_equal(obj3$logL, -1308.57082500809)

#expect_equivalent_to_reference(obj3, "spl2Dcylinder", tolerance = 1e-6)

## cyclic: torus
set.seed(1234)

sim_fun <- function(x) {
  y <- exp(sin(2*pi*x[1])) + exp(cos(2*pi*x[2]))
  y
}

# simulate data on a grid:
n1 <- 100
n2 <- 100
n <- n1*n2
dat <- expand.grid(x1=seq(0,1,length=n1),x2=seq(0,1,length=n2))
dat$ytrue <- apply(dat, MARGIN=1, FUN = sim_fun)

Ntraining <- 2500
k <- sample(x=c(1:n), size = Ntraining)
dat_train <- dat[k, ]

dat_train$y <- dat_train$ytrue + rnorm(n=Ntraining)

nseg <- c(20, 20)

obj4 <- LMMsolve(fixed = y~1,
                spline = ~spl2D(x1, x2, nseg, cyclic=c(TRUE, TRUE)),
                data = dat_train)

expect_equal(obj4$logL, -1287.57534355148)

#expect_equivalent_to_reference(obj4, "spl2Dtorus", tolerance = 1e-6)
