## Patch release

- Fixed check issues on macM1
- Minor improvements ready for release

----

## Test environments

* local Windows 10 install, R 4.1.3
* winbuilder (develop)
* macbuilder (release)
* macM1 builder
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Bart-Jan van Rossum <bart-jan.vanrossum@wur.nl>'

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

