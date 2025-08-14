f2 <- function(x) { 0.3 + 0.4*x + 0.2*sin(20*x) }

set.seed(12)
n <- 150
x <- seq(0, 1, length = n)
sigma2e <- 0.04
y <- f2(x) + rnorm(n, sd = sqrt(sigma2e))
dat <- data.frame(x, y)

obj <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 data = dat)

## Test for newdata
expect_error(predict(obj, newdata = "test", se.fit = TRUE),
             "newdata should be a data.frame")
newdat <- data.frame(x = seq(0, 1, length = 300))
newdat2 <- newdat
newdat2[1,] <- NA
expect_error(predict(obj, newdata = newdat2, se.fit = TRUE),
             "Variable x has missing values", fixed=TRUE)

newdat <- data.frame(x = seq(0, 1, length = 300))
pred1 <- predict(obj, newdata = newdat, se.fit = FALSE, type = "response")
pred2 <- predict(obj, newdata = newdat, se.fit = TRUE, type = "response")
pred3 <- predict(obj, newdata = newdat, se.fit = FALSE, type = "link")
pred4 <- predict(obj, newdata = newdat, se.fit = TRUE, type = "link")
pred5 <- predict(obj, newdata = newdat, se.fit = TRUE, deriv = "x")

expect_equivalent_to_reference(pred1, "pred1")
expect_equivalent_to_reference(pred2, "pred2")
expect_equivalent_to_reference(pred3, "pred3")
expect_equivalent_to_reference(pred4, "pred4")
expect_equivalent_to_reference(pred4, "pred5")

expect_error(predict(obj, newdata = newdat, se.fit = TRUE, deriv = "z2"),
             "Cannot find derivative for z2", fixed=TRUE)
expect_error(predict(obj, newdata = newdat, se.fit = TRUE, deriv = c("a","b")),
             "Argument deriv should have length one.")
expect_error(predict(obj, newdata = newdat, se.fit = TRUE, deriv=1),
             "Argument deriv should be a character string")


# multinomial example
set.seed(1234)

k <- 4
mu <- c(0.1, 0.4, 0.6, 0.9)
names(mu) <- LETTERS[1:k]

nonlinear <- function(x, mu) {
  z <- sapply(mu, function(mu) { exp(-8*sin(pi*(x-mu))^2)})
  # normalize to sum equal to one
  z <- z/sum(z)
  return(z)
}

x <- runif(n, 0, 1)
sz <- 10
multiNom <- t(sapply(x, FUN=
                       function(x) {
                         rmultinom(n=1, size=sz, prob = nonlinear(x,mu))
                       } ))
colnames(multiNom) <- names(mu)
dat <- data.frame(x, multiNom)
head(dat, 4)

obj <- LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
                spline = ~spl1D(x, nseg = 17, xlim = c(0,1)),
                data = dat,
                family = multinomial())
summary(obj)
sRows <- rowSums(multiNom)
fr <- multiNom/sRows
dat_fr <- data.frame(x, fr)

x0 <- seq(0, 1, by = 0.01)
newdat <- data.frame(x = x0)
expect_error(predict(obj, newdata = newdat, se.fit=TRUE),
             "se.fit=TRUE not implemented yet for multinomial")

pred5 <- predict(obj, newdata = newdat, se.fit=FALSE)
expect_equivalent_to_reference(pred5, "pred6")
