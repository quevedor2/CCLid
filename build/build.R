#install.packages("devtools")
library(devtools)

pkg=file.path("~/git", "CCL_authenticator")

#### Assembling data ####
usethis::use_data_raw()
devtools::load_all()

#### Building ####
setwd(pkg)
# create("CCLid")
devtools::document(pkg)
devtools::check(pkg)

#devtools::use_vignette("CCL-vignette", pkg)
#devtools::build_vignettes(pkg)
devtools::build(pkg)

#devtools::install_github("quevedor2/CCLid")
devtools::install("~/git/CCL_authenticator")
devtools::install("~/git/CCLid")
#devtools::reload(pkgload::inst('CCLid'))

devtools::install_github("quevedor2/CCL_authenticator")
