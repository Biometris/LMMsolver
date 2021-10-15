multipop <- read.csv("./data_raw/multipopQTL.csv")

multipop[["cross"]] <- as.factor(multipop[["cross"]])
multipop[["ind"]] <- as.factor(multipop[["ind"]])

colnames(multipop)[3:5] <- c("pA", "pB", "pC")

usethis::use_data(multipop, overwrite = TRUE)
