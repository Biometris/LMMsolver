#' Linear Mixed Model Solver using sparse matrix algebra.
#'
#' An efficient and flexible system to solve sparse mixed model
#' equations, for models that are often used in statistical genetics.
#' Important applications are the use of splines to model spatial or temporal
#' trends. Another application area is mixed model QTL analysis for
#' multiparental populations, allowing for heterogeneous residual variance and
#' random design matrices with Identity-By-Descent (IBD) probabilities.
#'
#' @name LMMsolver
#' @aliases LMMsolver LMMsolver-package
#' @docType package
#' @title Package LMMsolver
#' @author Martin Boer \email{martin.boer@@wur.nl}
#' @author Bart-Jan van Rossum \email{bart-jan.vanrossum@@wur.nl} (maintainer)
#' @references Martin P. Boer (2023). \emph{Tensor product P-splines using a sparse mixed model formulation},
#' Statistical Modelling, 23, p. 465 - 479. \url{https://doi.org/10.1177/1471082X231178591}
#'
#' @useDynLib LMMsolver, .registration = TRUE
#' @importFrom Rcpp sourceCpp
# The next and last line should be the word 'NULL'.
NULL
