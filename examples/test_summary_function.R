library(agridat)

# Test 1: simple standard mixed model

## Fit models on john.alpha data from agridat package.
data(john.alpha, package = "agridat")

## Fit the same model with genotype as fixed effect.
obj <- LMMsolve(fixed = yield ~ rep + gen, data = john.alpha)
summary(obj)
summary(obj, which='variances')

# baseline model, as in JABES 2020 paper, Table: 69.91:
round(deviance(obj), 2)

## Fit the same model with genotype as random effect.
obj <- LMMsolve(fixed = yield ~ rep,
                      random = ~gen,
                      data = john.alpha)
summary(obj)
summary(obj, which='variances')

N <- nrow(john.alpha)

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield~rep+gen,
                 spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord = 1, scaleX=FALSE),
                 data = john.alpha,
                 trace = FALSE,
                 tolerance = 1.0e-10)
summary(obj1)
# see Table 1:
round(deviance(obj1),2)

# Test 2: heterogeneous residual error:

## The residual variances for the two populations can be different.
## Allow for heterogeneous residual variances using the residual argument.
lGrp <- list(QTL = 3:5)
obj <- LMMsolve(fixed = pheno ~ cross,
                       group = lGrp,
                       random = ~grp(QTL),
                       residual = ~cross,
                       data = multipop,trace=TRUE)
summary(obj)
summary(obj, which='variances')

# Test 3: spl2D for SpATS SAP model:

# use durban data, as in BioRxiv 2021 paper Hans-Peter:
data(durban.rowcol)

dat <- durban.rowcol
dat <- dplyr::rename(dat, col=bed)
head(dat)

# Create factor variable for row and columns
dat$R <- as.factor(dat$row)
dat$C <- as.factor(dat$col)

obj <- LMMsolve(yield~rep,
               random=~R+C+gen,
               spline=~spl2D(col, row, nseg = c(20, 15)),
               data = dat)

summary(obj)
summary(obj, which='variances')
round(deviance(obj), 2)


