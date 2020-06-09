#install.packages("devtools")
library(devtools)

pkg=file.path("~/git", "CCLid")
pkg=file.path("~/git", "CCL_authenticator")

#### Assembling data ####
usethis::use_data_raw()
devtools::load_all()

#### Building ####
setwd(pkg)
Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")
# create("CCLid")
devtools::document(pkg)
devtools::check(pkg)


#usethis::use_vignette("CCL-vignette", pkg)
devtools::build_vignettes(pkg)

devtools::build(pkg)
#devtools::install_github("quevedor2/CCLid")
devtools::install("~/git/CCL_authenticator")
devtools::install("~/git/CCLid")
devtools::install(pkg)
#devtools::reload(pkgload::inst('CCLid'))

devtools::install_github("quevedor2/CCL_authenticator")
