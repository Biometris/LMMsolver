# APSIM analysis using mixed model P-splines
#
# Daniela Bustos-Korts and Martin Boer, Biometris, WUR.

## use
library(dplyr)

###############   Import data and specify input  ###############
data.HTP <- read.csv("./data_raw/HTP_BiomassOverTime.csv",
                     quote = "\'", stringsAsFactors = FALSE, na.string = "*")

# select the trait HERE!
pheno <- "biomass_MediumError"

# select genotype here
sel_geno <- "g001"

# define days between two measurements
freq <- 1

#################### some formatting #############################
# remove from beginning, only take observations after 20 days:
data.HTP <- filter(data.HTP, da_sow2 >= 20 )

# put in correct order
data.HTP <- arrange(data.HTP, geno, da_sow2)

# select measurement days
selected.days <- seq(min(data.HTP$da_sow2), (max(data.HTP$da_sow2)-1), freq)
data.HTP2 <-  as.data.frame(filter(data.HTP, da_sow2 %in% selected.days))

# Select genotype
data.HTP2_geno <- subset(data.HTP2, geno == sel_geno)

# some renaming.
dat <- data.HTP2_geno
dat <- dat %>% rename(das=da_sow2, biomass = biomass_MediumError) %>%
    select(env, geno, das, biomass)
APSIMdat <- dat

usethis::use_data(APSIMdat, overwrite = TRUE)
