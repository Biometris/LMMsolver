## Check with winbuilder - develop only.
devtools::check_win_devel()

## Check with macbuilder.
devtools::check_mac_release()

rhub::platforms()

## Check with rhub.
rhub::check_for_cran(path = "C:/Projects/R_packages/LMMsolver/")

## Check with rhub.
## Some specific platforms that give issues in earlier releases.
rhub::check_for_cran(platforms = c("ubuntu-gcc-release",
                                   "fedora-clang-devel",
                                   "debian-clang-devel",
                                   "debian-gcc-devel",
                                   "macos-highsierra-release"),
                     path = "C:/Projects/R_packages/LMMsolver/")


## For MAC M1 builder use this:
## This gave test issues for version 1.0.1
##https://mac.r-project.org/macbuilder/submit.html

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




