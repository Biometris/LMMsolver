## New release

- Initial CRAN release
- After earlier release fixed links in vignette and README

----

## Test environments

* local Windows 10 install, R 4.1.1
* winbuilder (release)
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Bart-Jan van Rossum <bart-jan.vanrossum@wur.nl>'

New submission

Possibly misspelled words in DESCRIPTION:
  IBD (9:54)
  multiparental (8:5)
  
These are spelled correctly.   

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.2307/1390762
    From: inst/doc/Solving_Linear_Mixed_Models.html
    Status: 403
    Message: Forbidden
  URL: https://www.jstor.org/stable/2246049
    From: inst/doc/Solving_Linear_Mixed_Models.html
    Status: 403
    Message: Forbidden

These links work correctly when opened from my browser.

