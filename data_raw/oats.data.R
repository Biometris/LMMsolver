library(agridat)
data("john.alpha", package = "agridat")
write.csv(john.alpha, file="oats.data.csv", quote=FALSE, row.names = FALSE)

oats.data <- read.csv("oats.data.csv", stringsAsFactors = TRUE)
usethis::use_data(oats.data, overwrite = TRUE)
