# Martin Boer, Biometris, Wageningen

# example sparse 1D, see Boer et al 2020, JABES for details.
rm(list=ls())

library(agridat)
library(LMMsolver)

data(john.alpha)
dat <- john.alpha
str(dat)

obj <- LMMsolve(fixed = yield~rep+gen,
                spatial = ~spl1D(x = plot, nseg = 20),
                data = dat,
                trace = TRUE,
                tolerance = 1.0e-8)
obj$logL
obj$ED

