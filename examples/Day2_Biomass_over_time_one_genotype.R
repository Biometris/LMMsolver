# APSIM analysis using mixed model P-splines
#
# Daniela Bustos-Korts and Martin Boer, Biometris, WUR.

## use
rm(list = ls())
library(dplyr)
#source("PsplinesREML.R")

###############   Import data and specify input  ###############
data.HTP<-read.csv('HTP_BiomassOverTime.csv',
         quote="\'",stringsAsFactors=FALSE, na.string="*")

# parameters to set......

# select the trait HERE!
pheno <- 'biomass_MediumError'

# select genotype here
sel_geno <-'g001'

# define days between two measurements
freq <- 1

#################### some formatting #############################
# remove from beginning, only take observations after 20 days:
data.HTP = filter(data.HTP,da_sow2 >= 20 )

# put in correct order
data.HTP = arrange(data.HTP,geno,da_sow2)

#select measurement days
selected.days =seq(min(data.HTP$da_sow2), (max(data.HTP$da_sow2)-1), freq)
data.HTP2 = as.data.frame(filter(data.HTP, da_sow2 %in% selected.days))

#Select genotype
data.HTP2_geno<- subset(data.HTP2,geno==sel_geno)

# some renameing.
dat <- data.HTP2_geno
dat <- dat %>% rename(das=da_sow2, biomass=biomass_MediumError) %>%
    select(env, geno, das, biomass)
head(dat)
nrow(dat)
APSIMdat <- dat

## save to file..
save(APSIMdat, file="APSIMdat.rda")

# run
obj2 <- LMMsolve(biomass~1,spline=~spl1D(x=das,nseg=200),
                 data=dat)
summary(obj2)

