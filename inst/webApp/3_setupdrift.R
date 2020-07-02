################################
#### Section 3: Drift setup ####
# 3_setupdrift.R
library(CCLid)
args = commandArgs(trailingOnly=TRUE)
tmp.dir <- args[1]
tmp.id <- args[2]

load(file.path(tmp.dir, paste0("tmp2_pred.", tmp.id, ".rda")))  #Starts with 82Mb RES
metadata <- c("meta.df",  "snp6.dat")
x <- sapply(metadata, downloadRefCCL, saveDir=sys.data, verbose=verbose)

## Look for drift
all.ids <- unique(unlist(pred[,c('Var1', 'Var2')]))
ref.dat <- CCLid::loadRef(sys.data, 'baf', bin.size=bin.size, 
                          meta.df=meta.df, just.var=TRUE)
vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, 
                      ref.mat=ref.dat$ref, ids=all.ids,
                      snp6.dat=snp6.dat)  ## starts at 2.5g, ramps up to 6.8Gb
save(bin.size, sys.data, vcfFile, vcf.mat, 
     file=file.path(tmp.dir, paste0("tmp3_drift_vcf.", tmp.id, ".rda")))


