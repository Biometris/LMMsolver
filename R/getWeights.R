getWeights <- function(weights, data) {
  if (!is.null(weights)) {
    if (!(weights %in% colnames(data))) {
      stop("weights not defined in dataframe data")
    }
    w <- data[[weights]]
    if (!is.numeric(w) || sum(is.na(w)) != 0 || min(w) < 0) {
      stop("weights should be a numeric vector with non-negative values")
    }
  } else {
    w <- rep(1, nrow(data))
  }
  return(w)
}
