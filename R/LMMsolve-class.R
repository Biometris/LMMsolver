#' Class LMMsolve
#'
#' Class LMMsolve of fitted mixed models.
#'
#' @section Slots:
#' An object of class LMMsolve contains the following slots.
#' \describe{
#' \item{logL}{The restricted log-likelihood at convergence}
#' \item{sigma2e}{The residual error}
#' \item{tau2e}{The estimated variance components}
#' \item{ED}{The effective dimensions}
#' \item{EDmax}{The maximal effective dimensions}
#' \item{EDnames}{The names of the effective dimensions}
#' \item{a}{The estimated effects from the mixed model equations}
#' \item{yhat}{The fitted values}
#' \item{dim}{The dimensions for each of the fixed and random terms in the
#' mixed model}
#' \item{term.labels}{The Names of the fixed and random terms in the mixed
#' model}
#' \item{splRes}{An object with definition of spline argument}
#' }
#'
#' @usage NULL
#'
#' @rdname LMMsolve-class
#' @export
new_LMMsolve <- function(object) {
  structure(object,
            class = c("LMMsolve", "list"))
}
