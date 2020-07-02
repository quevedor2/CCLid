##################################
#### Section 5: WebApp Return ####
# 5_return.R
library(CCLid)
args = commandArgs(trailingOnly=TRUE)
tmp.dir <- args[1]
tmp.id <- args[2]

load(file.path(tmp.dir, paste0("tmp4_drift.", tmp.id, ".rda")))  #Starts with 72Mb RES
load(file.path(tmp.dir, paste0("tmp2_pred.", tmp.id, ".rda")))  #Starts with 72Mb RES

## Cleanup:
genfiles <- list.files(tmp.dir)
if(all(c(paste0("seg.", tmp.id, ".json"), 
         paste0("drift.", tmp.id, ".json"), 
         paste0("pred.", tmp.id, ".json")) %in% genfiles)){
  setwd(tmp.dir)
  file.remove(c(paste0('tmp1_vcf_all.', tmp.id, '.rda'), 
                paste0('tmp2_pred.', tmp.id, '.rda'),
                paste0('tmp3_drift_vcf.', tmp.id, '.rda'), 
                paste0('tmp4_drift.', tmp.id, '.rda')))
}

dat <- list("pred"=pred,
            "seg"=bdf$cna.obj[[1]]$output,
            "fraction"=bdf$frac[[1]])
return(dat)
