#' Function to get the Effective Dimensions.
#'
#' @param object an object of class LMMsolve
#' @returns A data.frame with the effective dimensions and penalties.
#'
#'#' @examples
#' ## Fit model on oats data
#' data(oats.data)
#'
#' ## Fit a model with a 1-dimensional spline at the plot level.
#' obj <- LMMsolve(fixed  = yield ~ rep + gen,
#'                 spline = ~spl1D(x = plot, nseg = 20),
#'                 data = oats.data)
#' effDim(obj)
#
#' @export
effDim <- function(object) {
  return(object$EDdf)
}
