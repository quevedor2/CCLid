###################################
#### Section 1: Identity setup ####
# 1_setupid.R
library(CCLid)
args = commandArgs(trailingOnly=TRUE)
tmp.dir <- args[1]
vcfFile <- args[2]
bin.size <- args[3]
sys.data <- args[4]
num.snps <- args[5]
tmp.id <- args[6]

dir.create(tmp.dir, recursive=TRUE, showWarnings = FALSE)
metadata <- c("meta.df",  "snp6.dat")
x <- sapply(metadata, downloadRefCCL, saveDir=sys.data, verbose=verbose)

ref.dat <- CCLid::loadRef(sys.data, 'baf', bin.size=bin.size, 
                          just.var=TRUE, meta.df=meta.df)

## Load in VCF file of external data
vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, 
                      ref.mat=ref.dat$ref, max.snps=as.integer(num.snps),
                      snp6.dat=snp6.dat)
save(vcf.mat, bin.size, sys.data, vcfFile, 
     file=file.path(tmp.dir, paste0("tmp1_vcf_all.", tmp.id, ".rda"))) ## Ends with 92Mb RES

