# 
# option_list <- list(
#   make_option(c("-r", "--refdir"), type="character", default='CCLE',
#               help="Dataset to use, either 'GDSC' or 'GNE' or 'CCLE'"),
#   make_option(c("-v", "--vcf"), type="integer", default=10,
#               help="Size of the groups to run in one job [default= %default]", metavar="integer"))
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)
# 
# 
# PDIR <- opt$refdir
# vcfFile <- opt$vcf
# 
# ## Load in Ref mat file
# ref.dat <- CCLid::loadRef(PDIR, 'baf', bin.size=5e5)
# 
# ## Load in VCF file of external data
# vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, ref.mat=ref.dat$ref)
# 
# ## Look for similarity
# sample=basename(vcfFile)
# colnames(vcf.mat)[1] <- sample
# pred <- checkForConcordance(vcf.mat, sampleID=sample) 
# 
# ## Look for drift
# all.ids <- unique(unlist(pred$pred$M[,c('Var1', 'Var2')]))
# bdf <- bafDrift(vcf.mat[,c(sample, all.ids[grep(sample, all.ids, invert = TRUE)])])
# #CCLid:::plot.CCLid(bdf$cna.obj[[2]])
# 
# ## Return finished object for WebApp
# list("pred"=pred$pred$M,
#      "drift"=bdf$frac[[1]],
#      "seg"=bdf$cna.obj[[sample]]$output)
