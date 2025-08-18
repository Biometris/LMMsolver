set.seed(1234)

n <- 100
nobs <- 1000

fun_lambda <- function(x) { 0.5 + 0.4*sin(2*pi*x) }

x <- seq(0, 1, length=n)
dat <- data.frame(nr=sample(1:n, size = nobs, replace=TRUE))
time.dict <- data.frame(nr=1:n, x=x)
dat <- merge(dat, time.dict, by='nr', all.x=TRUE)
dat$score <- sapply(dat$x, FUN = function(x) {
                                    rbinom(n = 1, size = 1, prob = fun_lambda(x))
                                 })
dat$scoreF <- as.factor(ifelse(dat$score == 0,"failure","succes"))

# Option 1 for response, see generic ?family

mod1 <- LMMsolve(fixed = scoreF ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 data = dat)

tab <- table(dat$nr, dat$score)
dat2 <- data.frame(succes = tab[,2], failure=tab[,1])
dat2$cnt <- dat2$failure + dat2$succes
dat2$x <- x
dat2$y <- dat2$succes/dat2$cnt

# Option 2 for response, see generic ?family

mod2 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 weights = "cnt",
                 data = dat2)

# Option 3 for response, see generic ?family

mod3 <- LMMsolve(fixed = cbind(succes, failure) ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 data = dat2)

# MB, Aug 18 2025
# set family to NULL, change in r-devel of binomial function:
# okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")  <- 4.5.1
# okLinks <- c("logit", "probit", "cloglog", "cauchit", "log", "identity") <- new
#
mod1$family <- NULL
mod2$family <- NULL
mod3$family <- NULL
expect_equivalent_to_reference(mod1, "binomial1")
expect_equivalent_to_reference(mod2, "binomial2")
expect_equivalent_to_reference(mod3, "binomial3")

expect_error(LMMsolve(fixed = cbind(succes, failure, cnt) ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 data = dat2),
             "family binomial : response should have two columns.")

expect_error(LMMsolve(fixed = cbind(succes, failure) ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 weights = "cnt",
                 data = dat2),
      "family binomial : weights cannot be used in cbind format.")

dat3 <- dat2
dat3$succes[1] <- NA

expect_error(LMMsolve(fixed = cbind(succes, failure) ~ 1,
                      spline = ~spl1D(x, nseg = 50),
                      family = binomial(),
                      data = dat3),
      "family binomial : NA's in response variable.")

dat4 <- dat2
dat4$succes[1] <- "A"
expect_error(LMMsolve(fixed = cbind(succes, failure) ~ 1,
                      spline = ~spl1D(x, nseg = 50),
                      family = binomial(),
                      data = dat4),
      "response cbind(succes, failure) should be numeric.", fixed = TRUE)

dat5 <- dat2
dat5$succes[1] <- -1
expect_error(LMMsolve(fixed = cbind(succes, failure) ~ 1,
                      spline = ~spl1D(x, nseg = 50),
                      family = binomial(),
                      data = dat5),
             "family binomial : negative values in response variable.")






