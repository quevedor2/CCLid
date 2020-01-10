#install.packages("devtools")
library(devtools)

pkg=file.path("~/git", "CCL_authenticator")

devtools::document(pkg)
devtools::check(pkg)

#devtools::use_vignette("CCL-vignette", pkg)
#devtools::build_vignettes(pkg)
devtools::build(pkg)

devtools::install_github("quevedor2/CCL_authenticator")
