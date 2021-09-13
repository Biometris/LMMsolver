#' ---
#' title: SAP example sparse using USprecip example
#' author: Martin Boer, Biometris, WUR, Wageningen, the Netherlands
#' ---

rm(list=ls())

library(spam)
library(dplyr)
library(SpATS)
library(LMMsolver)
library(ggplot2)
library(JOPS)
library(maps)
library(sp)
library(reshape2)

# Get USA outline
usa = map("usa", region = "main", plot = F)

# Get precipitation data from spam
data(USprecip)
dat = data.frame(USprecip)
# only use observed data
dat = subset(dat, infill==1)
nrow(dat) # 5906 true records, as in SAP2014 paper.

#' SAP model using SpATS.nogeno
#' ==================

nseg <- c(41, 41)
degr <- 3
pord <- 2
tolerance <- 1.0e-6

# Fit PS-ANOVA model
obj1 = SpATS.nogeno(response = "anomaly",
                   spatial = ~SAP(lat, lon,  nseg = nseg, pord=pord, degree=degr),
                   data = dat,
                   control = list(maxit = 100, tolerance = tolerance,
                                  traceing = 2, update.psi = FALSE))

# effective dimensions close to SAP2014 paper:
summary(obj1)

# Extract the surfaces
Trend = obtain.spatialtrend(obj1, grid = c(200, 200))

# Turn matrix into a "long" data frame
Mu = Trend$fit
rownames(Mu) = Trend$row.p
colnames(Mu) = Trend$col.p
dens = reshape2::melt(Mu)
names(dens) = c('x', 'y', 'z')

# Plot fit with contours
sl = T
ccol = 'blue'
capsize = 15

# Find points within US boundary for clipping
v = point.in.polygon(dens$x, dens$y, usa$x, usa$y)
dens = subset(dens, v == 1)

plt1 = ggplot(dens,  aes(x, y, fill = z)) +
  geom_raster(show.legend = sl) +
  scale_fill_gradientn(colours = terrain.colors(100))+
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Precipitation anomaly, smoothed and clipped') +
  geom_contour(aes(z = z), color = ccol, show.legend = T) +
  JOPS_theme() +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = capsize),
        axis.title.x = element_text(size = capsize),
        axis.title.y = element_text(size = capsize))
print(plt1)

#' solve using sparse formulation SAP
#' ==================

x1 <- dat$lat
x2 <- dat$lon
y  <- dat$anomaly

s <- proc.time()[3]

knots1 <- PsplinesKnots(min(x1)-1.0e-10, max(x1)+1.0e-10, degree=degr, nseg=nseg[1])
knots2 <- PsplinesKnots(min(x2)-1.0e-10, max(x2)+1.0e-10, degree=degr, nseg=nseg[2])

B1 <- as.spam(Bsplines(knots1, x1))
B2 <- as.spam(Bsplines(knots2, x2))
q1 <- ncol(B1)
q2 <- ncol(B2)

D1 <- diff.spam(diag(q1), diff=pord)
DtD1 <- crossprod(D1)

D2 <- diff.spam(diag(q2), diff=pord)
DtD2 <- crossprod(D2)

# we have to calculate RowKronecker product only once:
B12 <- as.spam(RowKronecker(B1, B2))

# null spaces of DtD1 and DtD2 (not normalized!)
if (pord == 1)
{
  U1_null <- matrix(data=1, ncol=pord, nrow=q1)
  U2_null <- matrix(data=1, ncol=pord, nrow=q2)
} else {
  # L2 norm of a vector:
  norm_vec <- function(x) { return(sqrt(sum(x^2)))}

  # calculate the linear/fixed parts:
  U1_null <- cbind(1, scale(1:q1))
  U2_null <- cbind(1, scale(1:q2))

  U1_null <- apply(U1_null, MARGIN=2, function(x) (x/norm_vec(x)))
  U2_null <- apply(U2_null, MARGIN=2, function(x) (x/norm_vec(x)))
}

U_null <- kronecker(U1_null, U2_null)

X <- B12 %*% U_null
# take first column as intercept...
X[,1] <- 1

# make Rinv
lRinv = list()
Nobs <- length(y)
lRinv[[1]] = diag.spam(1, Nobs)
names(lRinv) = "residual"

# sparse boundary constraints
if (pord == 1)
{
  # other constraints also possible....
  C <- spam(x=c(1, rep(0,q1*q2-2), 1), nrow=q1*q2, ncol=1)
} else {
  C1 <- spam(x=0, nrow=q1, ncol=pord)
  C1[1,1] = C1[q1,2] = 1
  C2 <- spam(x=0, nrow=q2, ncol=pord)
  C2[1,1] = C2[q2,2] = 1
  C <- kronecker(C1, C2)
}

CCt <- as.spam(C %*% t(C))
lGinv <- list()
lGinv[[1]] <- kronecker(DtD1, diag.spam(q2)) + CCt
lGinv[[2]] <- kronecker(diag.spam(q1), DtD2) + CCt
names(lGinv) <- c('f(lat,lon)|lat', 'f(lat,lon)|lon')

obj2 = sparseMixedModels(y, X, B12, lGinv, lRinv,
              maxit=100, tolerance=tolerance, display=TRUE, trace=TRUE)
e <- proc.time()[3]
cat("Computation time:", e-s, "seconds")
obj2$ED

# compare efffective dimensions...
obj1$eff.dim[c(5,6)] - obj2$ED[c(2,3)]

t(C) %*% obj2$a[-c(1:pord^2)]

#' make predictions on a grid and compare
#' ==================

x1grid <- seq(min(x1), max(x1), length=200)
x2grid <- seq(min(x2), max(x2), length=200)
Bx1 <- as.spam(Bsplines(knots1, x1grid))
Bx2 <- as.spam(Bsplines(knots2, x2grid))
B12x <- Bx1 %x% Bx2
Xpred <- B12x %*% U_null
Xpred[,1] <- 1.0

# fit obtained from SpATS.nogeno,
# for comparison include intercept:
mu <- obj1$coeff[1]
fit1 = Trend$fit + mu

a <- obj2$a
bc <- as.vector(Xpred %*% a[1:pord^2])
sc <- as.vector(B12x %*% a[-c(1:pord^2)])
fit2 <- bc + sc
range(fit1 - fit2)

#' compare deviances
#' ==================

# compare SAP SpATS with LMMsolver..
dev1 <- as.numeric(obj1$deviance)
dev2 <- -2*obj2$logL
dev1
dev2

dev1-dev2

