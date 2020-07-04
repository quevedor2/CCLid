## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install-pkg, eval=FALSE, results='hide'----------------------------------
#  library(devtools)
#  devtools::install_github('bhklab/CCLid')

## ---- eval=FALSE, results='hide'----------------------------------------------
#  library(CCLid)

## ----download-dataset, eval=FALSE, results='hide'-----------------------------
#  refdir <- '~/cclid'
#  metadata <- c("meta.df",  "affy.omni", "cin70", "gne.meta", "melt.cells", "snp6.dat")
#  sapply(metadata, downloadRefCCL, saveDir=refdir, verbose=verbose)
#  
#  ref_dat <- CCLid::loadRef(PDIR=refdir, analysis='baf',
#                            bin.size=5e5, just.var=TRUE)

## ----load-vcf, eval=FALSE, results='hide'-------------------------------------
#  path_to_vcf = file.path(system.file(file.path("extdata"), package="CCLid"), "A549_trim.vcf.gz")
#  vcf_map <- mapVcf2Affy(path_to_vcf)

## ----map-vcf, eval=FALSE, results='hide'--------------------------------------
#  path_to_vcf = file.path(system.file(file.path("extdata"), package="CCLid"), "a549.sample_id.vcf")
#  vcf_mat <- compareVcf(path_to_vcf, var.dat=ref_dat$var, ref.mat=ref_dat$ref, max.snps=200, snp6.dat=snp6.dat)

## ----match-it, eval=FALSE, results='hide'-------------------------------------
#  colnames(vcf_mat)[1] <- sample_name
#  pred <- checkForConcordance(x.mat=vcf_mat, sampleID=sample_name, rm.gcsi = FALSE, meta.df=meta.df)

## ----match-it-`, eval=FALSE, results='hide'-----------------------------------
#  # Get all isogenic cell line IDs
#  all.ids <- unique(unlist(pred$pred$M[,c('Var1', 'Var2')]))

## ----match-it-2, eval=FALSE, results='hide'-----------------------------------
#  # Subset the VCF matrix for just the cell lines of interest
#  vcf_mat <- compareVcf(vcfFile, var.dat=ref.dat$var,
#                        ref.mat=ref.dat$ref, ids=all.ids,
#                        snp6.dat=snp6.dat)   ## starts at 2.5g, ramps up to 6.8Gb

## ----match-it-3, eval=FALSE, results='hide'-----------------------------------
#  # Drift estimation
#  bdf <- tryCatch({
#    bafDrift(vcf_mat, centering='median', snp6.dat=snp6.dat)
#  }, error=function(e){NULL})

