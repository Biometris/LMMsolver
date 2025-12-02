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

newdat2 <- data.frame(x= seq(-1,1, length = 10))
expect_error(predict(obj, newdata = newdat2, se.fit = TRUE, deriv="x"),
            "Variable x outside range of B-splines basis")

newdat3 <- data.frame(z = seq(-1,1, length = 10))
expect_error(predict(obj, newdata=newdat3),
             "variables (x) in data.frame newdata missing.", fixed=TRUE)

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

# make predictions for oats.data:
data(oats.data)
obj7 <- LMMsolve(yield ~ rep + gen,
                 random = ~block:rep,
                 data = oats.data)
newdat7 <- data.frame(rep = "R1", gen = levels(oats.data$gen))
pred7 <- predict(obj7, newdata = newdat7)
expect_equivalent_to_reference(pred7, "pred7")

# test if all fixed effects are defined in newdata:
newdat8 <- data.frame(gen = "G11")
expect_error(predict(obj7, newdata=newdat8),
                     "variables (rep) in data.frame newdata missing.", fixed=TRUE)

# example deriv for two dimensional splines

set.seed(1234)
f <- function(x,y) { x^2*y + 3*y^2*x}
df_dx <- function(x,y) { 2*x*y + 3*y^2}
df_dy <- function(x,y) { x^2 + 6*x*y }

n <- 1000
x <- runif(n, min=-1, max=1)
y <- runif(n, min=-1, max=1)
z <- f(x, y) + rnorm(n, sd = 0.01)
dat <- data.frame(x=x, y=y, z=z)

obj8 <- LMMsolve(fixed = z ~ 1,
                spline = ~spl2D(x1 = x, x2 = y, x1lim = c(-1,1), x2lim=c(-1,1),
                                nseg = c(25, 25)),
                data = dat)

newdat8 <- expand.grid(x = seq(-1, 1, length = 10),
                      y = seq(-1, 1, length = 10))

pred8 <- predict(obj8, newdata = newdat8, deriv = "x")
expect_equivalent_to_reference(pred8, "pred8")

pred9 <- predict(obj8, newdata = newdat8, deriv = "y")
expect_equivalent_to_reference(pred9, "pred9")

obj9 <- LMMsolve(fixed = z ~ 1,
                 spline = ~spl1D(x = x, nseg = 20,xlim=c(-1,1)) +
                           spl1D(x = y, nseg = 20,xlim=c(-1,1)),
                 data = dat)
expect_error(predict(obj9, newdat = newdat8, deriv = "x"),
             "Derivatives for multiple splines not implemented yet")


