library(ggplot2)
library(LMMsolver)

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
str(dat)
head(dat)
dim(dat)

obj0 <- LMMsolve(fixed = score ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 data = dat)
summary(obj0)

newdat <- data.frame(x = seq(0, 1, length = 300))
pred0 <- predict(obj0, newdata = newdat, se.fit=TRUE)
pred0$y_true <- fun_lambda(pred0$x)

# Option 1 for response, see generic ?family

obj1 <- LMMsolve(fixed = scoreF ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 data = dat)
summary(obj1)

pred1 <- predict(obj1, newdata = newdat, se.fit=TRUE)
pred1$y_true <- fun_lambda(pred1$x)

tab <- table(dat$nr, dat$score)
dat2 <- data.frame(succes = tab[,2], failure=tab[,1])
dat2$cnt <- dat2$failure + dat2$succes
dat2$x <- x
dat2$y <- dat2$succes/dat2$cnt

# Option 2 for response, see generic ?family

obj2 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 weights = "cnt",
                 data = dat2)
summary(obj2)

pred2 <- predict(obj2, newdata = newdat, se.fit = TRUE)
pred2$y_true <- fun_lambda(pred2$x)

# Option 3 for response, see generic ?family

obj3 <- LMMsolve(fixed = cbind(succes, failure) ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 data = dat2)
summary(obj3)

pred3 <- predict(obj3, newdata = newdat, se.fit=TRUE)
pred3$y_true <- fun_lambda(pred3$x)

ggplot(data = dat2, aes(x = x, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred3, aes(y = y_true), color = "red",
            linewidth = 1, linetype ="dashed") +
  geom_line(data = pred3, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data= pred3, aes(x=x, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  theme_bw()

# compare predictions.
all.equal(pred0, pred1)
all.equal(pred1, pred2)
all.equal(pred2, pred3)


