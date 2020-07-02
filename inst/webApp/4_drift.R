##########################
#### Section 4: Drift ####
# 4_drift.R
library(CCLid)
library(jsonlite)
args = commandArgs(trailingOnly=TRUE)
tmp.dir <- args[1]
outdir <- args[2]
tmp.id <- args[3]
verbose=FALSE

load(file.path(tmp.dir, paste0("tmp3_drift_vcf.", tmp.id, ".rda")))  #Starts with 72Mb RES
metadata <- c("snp6.dat")
x <- sapply(metadata, downloadRefCCL, saveDir=sys.data, verbose=verbose)

bdf <- bafDrift(vcf.mat, centering='median', snp6.dat=snp6.dat)
save(bdf, file=file.path(tmp.dir, paste0("tmp4_drift.", tmp.id, ".rda")))

if(verbose) cat("Writting to drift.json...")
sink(file=file.path(outdir, paste0("drift.", tmp.id, ".json")))
jsonlite::toJSON(lapply(bdf$frac[[1]], function(i){
  data.frame("frac"=i,
             "effect.size"=names(i))}), pretty = TRUE)
sink()

if(verbose) cat("Writting to seg.json...")
sink(file=file.path(outdir, paste0("seg.", tmp.id, ".json")))
jsonlite::toJSON(bdf$cna.obj[[1]]$output, pretty = TRUE)
sink()

## Cleanup:
# genfiles <- list.files(tmp.dir)
# if(all(c("seg.json", "drift.json", "pred.json") %in% genfiles)){
#   setwd(tmp.dir)
#   file.remove(c('tmp1_vcf_all.rda', 'tmp2_pred.rda', 'tmp3_drift_vcf.rda'))
# }




