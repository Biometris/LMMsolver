#' Simulated QTL mapping data set
#'
#' @format A data.frame with 180 rows and 6 columns.
#' \describe{
#' \item{cross}{Cross ID, two populations, AxB and AxC}
#' \item{ind}{Genotype ID}
#' \item{pA}{Probability that individual has alleles from parent A}
#' \item{pB}{Probability that individual has alleles from parent B}
#' \item{pC}{Probability that individual has alleles from parent C}
#' \item{pheno}{Simulated phenotypic value}
#' }
"multipop"

#' Simulated Biomass as function of time using APSIM wheat.
#'
#' @format A data.frame with 121 rows and 4 columns.
#' \describe{
#' \item{env}{Environment, Emerald in 1993}
#' \item{geno}{Simulated genotype g001}
#' \item{das}{Days after sowing}
#' \item{biomass}{Simulated biomass using APSIM; medium measurement error added}
#' }
#'
#' @references Bustos-Korts et al. (2019) Combining Crop Growth Modeling and
#' Statistical Genetic Modeling to Evaluate Phenotyping Strategies
#' \doi{10.3389/FPLS.2019.01491}
"APSIMdat"

#' Sea Surface Temperature
#'
#' @format A data.frame with 15607 rows and 4 columns.
#' \describe{
#' \item{lon}{longitude}
#' \item{lat}{latitude}
#' \item{sst}{sea surface temperature in Kelvin}
#' \item{type}{defines training and test set}
#' }
#'
#' @references Cressie et al. (2022) Basis-function models in spatial statistics.
#' Annual Review of Statistics and Its Application.
#' \doi{10.1146/annurev-statistics-040120-020733}
"SeaSurfaceTemp"

#' Uniformity trial of barley
#' @format A data.frame with 1076 rows and 3 columns
#' \describe{
#'  \item{row}{row coordinate}
#'  \item{col}{column coordinate}
#'  \item{yield}{yield per plot}
#' }
#' @source H. P. Piepho & E. R. Williams (2010). Linear variance models
#' for plant breeding trials. Plant Breeding, 129, 1-8.
#' \doi{/doi.org/10.1111/j.1439-0523.2009.01654.x}
#' @references Piepho, Hans‐Peter, Martin P. Boer, and Emlyn R. Williams.
#' "Two‐dimensional P‐spline smoothing for spatial analysis of plant breeding trials."
#' Biometrical Journal 64, no. 5 (2022): 835-857.
"barley.uniformity.trial"

#' Alpha lattice design of spring oats
#' @format A data.frame with 72 rows and 7 columns
#' \describe{
#' \item{plot}{plot number}
#' \item{rep}{replicate}
#' \item{block}{incomplete block}
#' \item{gen}{genotype}
#' \item{yield}{dry matter yield}
#' \item{row}{row}
#' \item{col}{column}
#' }
#' @details The response is grain yield in kg per hectare. The design was an alpha design
#' with 24 varieties, three replicates and  six incomplete blocks of size four per replicate.
#' The 72 plots were arranged in a single linear array.
#' @source J. A. John & E. R. Williams (1995). Cyclic and computer
#' generated designs. Chapman and Hall, London. Page 146.
#' @references Boer, Martin P., Hans-Peter Piepho, and Emlyn R. Williams.
#' "Linear variance, P-splines and neighbour differences for spatial adjustment
#' in field trials: how are they related?." JABES 25, no. 4 (2020): 676-698.
"oats.data"
