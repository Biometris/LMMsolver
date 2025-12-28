#' Heritability of genotype means
#'
#' Computes the heritability of genotype means (\eqn{h^2_{M_G}}) from a fitted
#' linear mixed model. This measure is equivalent to the SpATS, Cullis, and Oakey
#' heritability definitions used in modern mixed-model analyses of field trials.
#'
#' The heritability is defined as the effective dimension of the genotype term
#' divided by the maximum effective dimension of that term, assuming independent
#' genotypic effects. It reflects the reliability of predicted genotype means
#' rather than a variance ratio of raw phenotypes.
#'
#' @param obj An object of class \code{"LMMsolve"}.
#' @param geno.term A character string giving the name of the genotype term in
#'   the model. The term must be specified as random.
#'
#' @return A numeric value giving the heritability of genotype means
#'   (\eqn{h^2_{M_G}}).
#'
#' @details
#' This function returns a single heritability measure. Classical ANOVA-based
#' heritability estimates are not provided, as they are generally not meaningful
#' for unbalanced, incomplete block, or spatial designs.
#'
#' The returned value corresponds to
#' \deqn{h^2_{M_G} = \mathrm{ED}_g / \mathrm{ED}_{\max},}
#' where \eqn{\mathrm{ED}_g} is the effective dimension corresponding to the
#' genotypic effects and \eqn{\mathrm{ED}_{\max}} is the maximum effective
#' dimension of the genotype term, equal to the number of genotypes minus one.
#'
#' @references
#' Rodríguez-Álvarez, M. X., Boer, M. P., van Eeuwijk, F. A., & Eilers, P. H. C.
#' (2018).
#' Correcting for spatial heterogeneity in plant breeding experiments with
#' P-splines.
#' \emph{Spatial Statistics}, 23, 52–71.
#'
#' Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006).
#' Joint modeling of additive and non-additive genetic line effects in single
#' field trials.
#' \emph{Theoretical and Applied Genetics}, 113, 809–819.
#'
#' Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006).
#' On the design of early generation variety trials with correlated data.
#' \emph{Journal of Agricultural, Biological, and Environmental Statistics}, 11,
#' 381–393.
#'
#' @examples
#' \dontrun{
#' data(oats.data)
#'
#' obj <- LMMsolve(
#'   fixed  = yield ~ rep,
#'   random = ~ gen + rep:block,
#'   data   = oats.data
#' )
#'
#' getHeritability(obj, geno.term = "gen")
#' }
#'
#' @export
getHeritability <- function(obj, geno.term) {

  if (!inherits(obj, "LMMsolve")) {
    stop("obj must be an object of class 'LMMsolve'")
  }

  if (!is.character(geno.term) || length(geno.term) != 1L) {
    stop("geno.term must be a single character string")
  }

  EDdf <- effDim(obj)

  if (!(geno.term %in% EDdf$Term)) {
    stop(paste(geno.term, "not defined in the model"))
  }

  penalty <- subset(EDdf, Term == geno.term)$Penalty
  if (penalty < 1.0e-15) {
    stop(paste(geno.term, "should be in the random term"))
  }

  subset(EDdf, Term == geno.term)$Ratio
}
