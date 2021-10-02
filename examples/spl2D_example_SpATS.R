# Martin Boer, Biometris
#
# comparison between SpATS using SAP and implementation in
# LMMsolver.
#
library(SpATS)
library(LMMsolver)
library(agridat)
library(dplyr)

# use durban data, as in BioRxiv 2021 paper Hans-Peter:
data(durban.rowcol)

dat <- durban.rowcol
dat <- rename(dat, col=bed)
head(dat)

# Create factor variable for row and columns
dat$R <- as.factor(dat$row)
dat$C <- as.factor(dat$col)

# set parameters
tol = 1.0e-8
nseg = c(20, 15)
grid = c(80, 100)

m0 <- SpATS(response = "yield",
            spatial = ~ SAP(col, row, nseg = nseg),
            genotype = "gen",
            fixed = ~rep,
            random = ~ R + C,
            genotype.as.random = TRUE,
            data = dat,
            control = list(tolerance = tol))
# Brief summary
summary(m0)
plot(m0)
dev0 <- m0$deviance

# degree and pord not defined yet
m1 <- LMMsolve(yield~rep,
               random=~R+C+gen,
               spline=~spl2D(col, row, nseg = nseg),
               data = dat,
               trace=TRUE,
               tolerance = tol)
dev1 <- m1$dev

dev0
dev1
dev0 - dev1

M0 <- SpATS::obtain.spatialtrend(m0, grid=grid)$fit

df <- obtainSmoothTrend(m1, grid=grid)
M1 <- t(matrix(data = df$ypred, nrow=grid[1], ncol= grid[2],byrow=TRUE))

# only a constant difference...
range(M0 - M1)

# give summary
summary(m1)
