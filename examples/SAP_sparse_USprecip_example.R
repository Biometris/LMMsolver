#' ---
#' title: SAP example sparse
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
dat = subset(dat, infill==1)

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
                                  monitoring = 2, update.psi = FALSE))
summary(obj1)

# Extract the surfaces
Tr = obtain.spatialtrend(obj1, grid = c(200, 200))

# Turn matrix into a "long" data frame
Mu = Tr$fit
rownames(Mu) = Tr$row.p
colnames(Mu) = Tr$col.p
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
B21 <- as.spam(RowKronecker(B2, B1))

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

X <- B21 %*% (U2_null %x% U1_null)
# take first column as intercept...
X[,1] <- 1

# make Rinv
lRinv = list()
Nobs <- length(y)
lRinv[[1]] = diag.spam(1, Nobs)
names(lRinv) = "residual"

# boundary constraints
if (pord == 1)
{
  # other constraints also possible....
  C <- spam(x=c(1, rep(0,q1*q2-2), 1), nrow=q1*q2, ncol=1)
} else {
  C1 <- spam(x=0, nrow=q1, ncol=pord)
  C1[1,1] = C1[q1,2] = 1
  C2 <- spam(x=0, nrow=q2, ncol=pord)
  C2[1,1] = C2[q2,2] = 1
  C <- kronecker(C2, C1)
}

CCt <- as.spam(C %*% t(C))
lGinv <- list()
lGinv[[1]] <- kronecker(diag.spam(q2), DtD1) + CCt
lGinv[[2]] <- kronecker(DtD2, diag.spam(q1)) + CCt
names(lGinv) <- c('f(E,C)|E', 'f(E,C)|C')

obj2 = sparseMixedModels(y, X, B21, lGinv, lRinv,
              maxiter=100, eps=tolerance, display=TRUE, monitor=TRUE)
e <- proc.time()[3]
cat("Computation time:", e-s, "seconds")
obj2$ED

t(C) %*% obj2$a[-c(1:pord^2)]

#' compare results
#' ==================

# compare SAP SpATS with LMMsolver..
dev1 <- as.numeric(obj1$deviance)
dev2 <- -2*obj2$logL
dev1
dev2

dev1-dev2

