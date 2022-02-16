library(LMMsolver)
library(agridat)
library(asreml)
library(dplyr)

data(john.alpha)
dat <- john.alpha

# Full data.
obj1 <- asreml(fixed = yield~rep+gen+rep:block, data = dat)

# error, fixed part is not full rank....
obj2 <- LMMsolve(fixed = yield~rep+gen+rep:block,
                 data = dat,
                 trace = TRUE,
                 tolerance = 1.0e-10)

coef(obj1, list = TRUE)$`rep:block`
coef(obj2)$`rep:block`

# one block missing....
dat2 <- dat %>%
  filter(!(block=='B1' & rep=='R3')) %>%
  filter(!(block=='B2' & rep=='R1')) %>%
  filter(!(block=='B6' & rep=='R2')) %>%
  filter(!(block=='B4' & rep=='R3')) %>%
  filter(!(block=='B4' & rep=='R1'))

table(dat2$rep, dat2$block)

obj1a <- asreml(fixed = yield~rep+gen+rep:block,
                data = dat2)

tst <- lm(yield~rep+gen+rep:block,
          data = dat2)
coef(tst)

# error, fixed part is not full rank....
obj2a <- LMMsolve(fixed = yield~rep+gen+rep:block,
                  data = dat2,
                  trace = TRUE,
                  tolerance = 1.0e-10)

obj3a <- lm(yield~rep+gen+rep:block,
            data = dat2)

c1a <- coef(obj1a, list = TRUE)$`rep:block`[,1]
coef(obj1a, list = TRUE)$`(Intercept)`
coef(obj1a, list = TRUE)$rep

c2a <- coef(obj2a)$`rep:block`
coef(obj2a, list = TRUE)$`(Intercept)`
coef(obj2a)$rep

c3a <- coef(obj3a)

c1a[names(c2a)] - c2a



