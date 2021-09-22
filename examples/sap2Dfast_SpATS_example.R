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

# tolerance
tol = 1.0e-8

m0 <- SpATS(response = "yield",
            spatial = ~ SAP(col, row, nseg = c(10,20), degree = 3, pord = 2),
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
               spatial=~sap2D(col, row, knots = c(10,20), scaleX=TRUE),
               data = dat,
               trace=TRUE,
               tolerance = tol)
dev1 <- -2.0*m1$logL

dev0
dev1
dev0-dev1

