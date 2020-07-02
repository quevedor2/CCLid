#################################
#### Section 2: CCL Identity ####
# 2_cclid.R
library(CCLid)
library(jsonlite)
args = commandArgs(trailingOnly=TRUE)
tmp.dir <- args[1]
outdir <- args[2]
tmp.id <- args[3]
verbose=FALSE

load(file.path(tmp.dir, paste0("tmp1_vcf_all.", tmp.id, ".rda")))  #Starts with 92Mb RES
metadata <- c("meta.df")
x <- sapply(metadata, downloadRefCCL, saveDir=sys.data, verbose=verbose)

## Look for similarity
sample=gsub(".vcf.*", "", basename(vcfFile))
colnames(vcf.mat)[1] <- sample
pred <- checkForConcordance(x.mat=vcf.mat, sampleID=sample, 
                            rm.gcsi = FALSE, meta.df=meta.df)
pred <- pred$pred$M
save(pred, bin.size, sys.data, vcfFile, 
     file=file.path(tmp.dir, paste0("tmp2_pred.", tmp.id, ".rda")))

## Writting JSON and RDS files for WebApp
if(verbose) cat("Writting to pred.json...")
sink(file=file.path(outdir, paste0("pred.", tmp.id, ".json")))
jsonlite::toJSON(pred, pretty = TRUE)
sink()

