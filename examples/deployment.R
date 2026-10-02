## Check with winbuilder - develop only.
devtools::check_win_devel()

## Check with macbuilder.
devtools::check_mac_release()

rhub::platforms()

## Check with rhub.
rhub::check_for_cran(path = "C:/Projects/R_packages/LMMsolver/")

## Check with rhub.
## Check for specific platforms.
## A list will come up options.
rhub::rc_submit()

## Rebuild readme.
devtools::build_readme()

## Submit to CRAN
devtools::release()

## Build site for local check.
pkgdown::clean_site()
pkgdown::build_site()

## Code coverage - local.
detach("package:LMMsolver", unload = TRUE)
covr::gitlab(type = "none", code = "tinytest::run_test_dir(at_home = TRUE)")
library(LMMsolver)

## Check reverse dependencies. (takes about 45 minutes).
detach("package:LMMsolver", unload = TRUE)
revdepcheck::revdep_reset()
revdepcheck::revdep_check()




