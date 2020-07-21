#' ---
#' title: Hierarchical mixed model with sum to zero constraints.
#' author: Martin Boer, Biometris, WUR, Wageningen, The Netherlands
#' ---

rm(list=ls())
library(LMMsolver)
library(dplyr)

set.seed(1234)

#' Simulate some data
#' ---

# 5 genotypes
Ngeno <- 5

# 8 plants per geno
NplantsGeno <- 8

# 10 observations per plant
NobsPlants <- 10

# total number of plants.
totNplants = Ngeno*NplantsGeno

N <- totNplants*NobsPlants

mu <- 5.0
geno_eff <- rnorm(Ngeno, sd=1.0)
plant_eff <- rnorm(Ngeno*NplantsGeno, sd=2.0)
res_eff <- rnorm(N,sd=1.0)

x=1:NobsPlants
overalltrend <- 0.1*(x-5)^2
trend <- rep(overalltrend,totNplants)

y <- mu + rep(geno_eff,each=NplantsGeno*NobsPlants) +
      rep(plant_eff, each=NobsPlants) + trend + res_eff
y <- mu + rep(geno_eff,each=NplantsGeno*NobsPlants) +
  as.vector(kronecker(plant_eff, overalltrend)) + trend + res_eff

plantLabel <- paste0("P",formatC(1:totNplants,width=2,flag="0",format="d"))
genoLabel <- paste0("G", formatC(1:Ngeno,width=2,flag="0",format="d"))
plantID <- rep(plantLabel,each=NobsPlants)
genoID <- rep(genoLabel,each=NobsPlants*NplantsGeno)
obsnr <- rep(1:NobsPlants,times=totNplants)

df <- data.frame(geno=genoID, plant=plantID, timepoints = obsnr, trend=trend, y=y)
head(df)

# do some filtering for testing, not all genotypes have same number
# of plants.
#df <- filter(df, ! (plantID %in% c('P07','P28','P39')))
#y <- df$y

#' Genotype to Plant matrix Q and the projection matrix PQ:
#' ------

# construct the matrix Q:
tmp <- df %>% dplyr::select(geno,plant) %>% unique()
Q <- model.matrix(~geno-1, tmp)
rownames(Q) <- tmp$plant
colnames(Q) <- genoLabel
head(Q, 15)
# vector, with number of subjects per group
m_grp <- colSums(Q)
totNplants = sum(m_grp)
N <- nrow(df)

#'  projection matrix, orthogonal to Q:
PQ = diag(totNplants) - Q %*% solve(crossprod(Q)) %*% t(Q)
# Q and PQ are orthogonal....
range(PQ %*% Q)

#' Analysis of the data without splines.
#' ---

# what do we have:
#   k groups, total m subject, m_grp subjects per group, s timepoints per
#   subject (assumed to be balanced), N = m*s total number of observations:

k <- ncol(Q)
m <- sum(m_grp)  # totNplants
s <- length(unique(df$timepoints))
N <- m*s

# fixed effect
X = matrix(data=1, nrow=N, ncol=1)

# design matrix for plant:
U = kronecker(diag(m), rep(1,s))

# random effect for genotype
Z1 = U %*% Q
# random effect for plants, conditional on genotype:
Z2 = U %*% PQ
Z = cbind(Z1, Z2)
ncol(Z)
qr(Z)$rank

lZ <- list()
lZ[[1]] = Z1
lZ[[2]] = Z2
Z <- do.call("cbind", lZ)

df_ext = cbind(df, Z)

lM <- ndxMatrix(df, lZ, c("genotypes","plants"))

obj1 <- LMMsolve(fixed=y~1, randomMatrices=lM,data=df_ext, display=TRUE)
obj1$logL
obj1$ED

# built-in constraint, orthogonal to X:
sum(coef(obj1)$genotypes)

# constraints due to overlap in Z1 and Z2:
t(Q) %*% coef(obj1)$plants

# second formulation (assuming balanced data)

lZ <- list()
I = diag(NplantsGeno)
J = matrix(data=1,ncol=NplantsGeno,nrow=NplantsGeno)
K = I - (1/NplantsGeno)*J
Usc <- eigen(K)$vectors[, -NplantsGeno]
# for all genotypes same, balanced:
UscG <- kronecker(diag(1,Ngeno), Usc)
lZ[[1]] = kronecker(UscG,rep(1,NobsPlants))
Z <- do.call("cbind", lZ)

df_ext = cbind(df, Z)

lM <- ndxMatrix(df, lZ, c("plants"))

obj2 <- LMMsolve(fixed=y~1, random=~geno, randomMatrices=lM,data=df_ext, display=TRUE)
obj1$logL
obj2$logL

obj1$ED
obj2$ED

# built-in constraint, orthogonal to X:
sum(coef(obj2)$geno)

# plant effects are deviations from genotypic mean:
plant_eff <- UscG %*% coef(obj2)$plants
t(Q) %*% plant_eff

