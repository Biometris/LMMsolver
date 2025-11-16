#library(agridat)
#data("piepho.barley.uniformity", package = "agridat")
#data("barley.uniformity.trial")
# Remove NA to prevent spurious warnings.
#dat <- piepho.barley.uniformity[!is.na(piepho.barley.uniformity[["yield"]]), ]
#write.csv(dat, file="barley.uniformity.trial.csv", quote=FALSE, row.names = FALSE)

barley.uniformity.trial <- read.csv("barley.uniformity.trial.csv")
usethis::use_data(barley.uniformity.trial, overwrite = TRUE)

