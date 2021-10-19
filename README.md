
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LMMsolver

[![](https://www.r-pkg.org/badges/version/LMMsolver)](https://www.r-pkg.org/pkg/LMMsolver)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/LMMsolver)](https://www.r-pkg.org/pkg/LMMsolver)
[![R-CMD-check](https://github.com/Biometris/LMMsolver/workflows/R-CMD-check/badge.svg)](https://github.com/Biometris/LMMsolver/actions?workflow=R-CMD-check)
[![codecov](https://codecov.io/gh/Biometris/LMMsolver/branch/master/graph/badge.svg)](https://codecov.io/gh/Biometris/LMMsolver)

The aim of the `LMMsolver` package is to provide an efficient and
flexible system to solve sparse mixed model equations, for models that
are often used in statistical genetics. One important application is
using splines to model spatial Boer, Piepho, and Williams (2020) or
temporal (Bustos-Korts et al. 2019) trends for phenotypic observations.
Another application area is mixed model QTL analysis for multiparental
populations, allowing for heterogeneous residual variance and random
design matrices with Identity-By-Descent (IBD) probabilities (Li et al.
2021).

## Installation

-   Install from CRAN:

``` r
install.packages("LMMsolver")
```

-   Install latest development version from GitHub (requires
    [remotes](https://github.com/r-lib/remotes) package):

``` r
remotes::install_github("Biometris/LMMsolver", ref = "develop", dependencies = TRUE)
```

## Example

``` r
library(LMMsolver)
library(ggplot2)

## Get precipitation data from spam
data(USprecip, package = "spam")

## Only use observed data
USprecip <- as.data.frame(USprecip)
USprecip <- USprecip[USprecip$infill == 1, ]
# 5906 true records, as in 2014 paper.
nrow(USprecip) 
#> [1] 5906
```

The two-dimensional P-spline can be defined with the `spl2D()` function:

``` r
obj3 <- LMMsolve(fixed = anomaly ~ 1,
                 spline = ~spl2D(x1 = lon, x2 = lat, nseg = c(41, 41)),
                 data = USprecip)
```

Here the same number of segments is used as in M. Rodríguez-Álvarez and
Lee (2014). The summary function gives a table with the effective
dimensions and the penalty parameters:

``` r
summary(obj3)
#> Table with effective dimensions and penalties: 
#> 
#>          term Effective Model Nominal     Ratio     Penalty
#> 1 (Intercept)    1.0000     1       1 1.0000000  0.00000000
#> 2        splF    3.0000     3       3 1.0000000  0.00000000
#> 3      s(lon)  302.6004  1936    1932 0.1566255  0.09679905
#> 4      s(lat)  409.0868  1936    1932 0.2117426  0.36223681
#> 5    residual 5190.3128  5906    5906 0.8788203 13.52577271
#> 
#>  Total Effective Dimension: 5906
```

The spatial trend for the precipation can now be plotted on the map of
the USA.

``` r
plotDat <- obtainSmoothTrend(obj3, grid = c(200, 300), includeIntercept = TRUE)
usa = maps::map("usa", regions = "main", plot = FALSE)
v <- sp::point.in.polygon(plotDat$lon, plotDat$lat, usa$x, usa$y)
plotDat <- plotDat[v == 1, ]

ggplot(plotDat, aes(x = lon, y = lat, fill = ypred)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradientn(colors = topo.colors(100))+
  labs(title = "Precipitation (anomaly) US April 1948", x = "Longitude", y = "Latitude") +
  coord_fixed() +
  theme(panel.grid = element_blank())
```

<img src="man/figures/README-Plot_USprecip-1.png" width="70%" />

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Boer2020" class="csl-entry">

Boer, Martin P., Hans Peter Piepho, and Emlyn R. Williams. 2020. “<span
class="nocase">Linear Variance, P-splines and Neighbour Differences for
Spatial Adjustment in Field Trials: How are they Related?</span>” *J.
Agric. Biol. Environ. Stat.* 25 (4): 676–98.
<https://doi.org/10.1007/S13253-020-00412-4>.

</div>

<div id="ref-Bustos-Korts2019" class="csl-entry">

Bustos-Korts, Daniela, Martin P. Boer, Marcos Malosetti, Scott Chapman,
Karine Chenu, Bangyou Zheng, and Fred A. van Eeuwijk. 2019. “<span
class="nocase">Combining Crop Growth Modeling and Statistical Genetic
Modeling to Evaluate Phenotyping Strategies</span>.” *Front. Plant Sci.*
10 (November). <https://doi.org/10.3389/FPLS.2019.01491/FULL>.

</div>

<div id="ref-Li2021" class="csl-entry">

Li, Wenhao, Martin P. Boer, Chaozhi Zheng, Ronny V. L. Joosen, and Fred
A. van Eeuwijk. 2021. “<span class="nocase">An IBD-based mixed model
approach for QTL mapping in multiparental populations</span>.” *Theor.
Appl. Genet. 2021* 1 (August): 1–18.
<https://doi.org/10.1007/S00122-021-03919-7>.

</div>

<div id="ref-Rodriguez-Alvarez2018" class="csl-entry">

Rodríguez-Álvarez, María Xosé, Martin P. Boer, Fred A. van Eeuwijk, and
Paul H. C. Eilers. 2018. “<span class="nocase">Correcting for spatial
heterogeneity in plant breeding experiments with P-splines</span>.”
*Spat. Stat.* 23 (March): 52–71.
<https://doi.org/10.1016/J.SPASTA.2017.10.003>.

</div>

<div id="ref-Rodriguez-Alvarez2014" class="csl-entry">

Rodríguez-Álvarez, MX, and DJ Lee. 2014. “<span class="nocase">Fast
smoothing parameter separation in multidimensional generalized
P-splines: the SAP algorithm</span>.” *Stat. …*.
<http://link.springer.com/article/10.1007/s11222-014-9464-2>.

</div>

</div>
