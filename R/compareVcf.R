#' compareVcf
#' @description Checks a VCF file(s) against reference dataset
#' to look for similarity to any known cell lines
#'
#' @param vcfFile Path to VCF file to check against reference datasets
#' @param var.dat Variance data (list)
#' @param ref.mat Reference matrix (matrix)
#' @param max.snps Max number of SNPs to reduce 
#' @param ids IDs
#' @param sampletype Strictly for labelling purposes
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#' @param ... Extra param
#' @param all.ids IDs to subset to
#'
#' @return Matrix: Containing reference matrix subsetted to common SNPs as 
#' the input VCF, as well as a left-joined VCF data
#' @export
#'
compareVcf <- function(vcfFile, var.dat, ref.mat, 
                       max.snps=1e6, ids=NULL, sampletype='RNA', all.ids=NULL, 
                       snp6.dat, ...){
  vcf.map <- CCLid::mapVcf2Affy(vcfFile, snp6.dat=snp6.dat)
  vcf.map <- .filt(vcf.map, ...) ## Memory: up to 1.8Gb 
  
  ## Combine matrices and reduce features
  ## Find the overlap between the COMParator and the REFerence
  vcf.map.var <- CCLid::mapVariantFeat(vcf.map, var.dat)
  vcf.to.use <- vcf.map.var
  ov.idx <- CCLid::overlapPos(comp = vcf.to.use$BAF,
                              ref=ref.mat, mapping = 'probeset')
  if(nrow(ov.idx) > max.snps){
    ov.idx <- ov.idx[order(ov.idx$ref)[1:max.snps],]
  }
  rm(vcf.map, vcf.map.var); gc()  ## Cleanup
  
  ## BIG memory sink: shoots up to 7Gb
  if(is.null(ids)){
    refm <- ref.mat
  } else {
    colidx <- which(colnames(ref.mat) %in% all.ids)
    message("Isolating for: ", paste(colnames(ref.mat)[colidx], collapse=", "))
    refm <- ref.mat[,colidx,drop=FALSE]
  }
  x.mat <- cbind(vcf.to.use$BAF$BAF[ov.idx$comp], 
                 refm[ov.idx$ref,])
  colnames(x.mat)[1] <- paste0(sampletype, "_", gsub(".vcf.*", "", basename(vcfFile)))
  
  if(storage.mode(ref.mat[,1]) == 'integer'){
    x.mat[,-1] <- x.mat[,-1] / 100
  }
  return(x.mat)
}

