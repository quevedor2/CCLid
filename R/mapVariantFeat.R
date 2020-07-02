#' mapVariantFeat
#' @description Selects only the features with the highest variance
#' per bin
#' 
#' @param vcf.map An object returned by mapVcf2Aff()
#' @param var.dat A list returned by formatRefMat()$var
#'
#' @export
mapVariantFeat <- function(vcf.map, var.dat){
  is.bm <- FALSE
  if(class(vcf.map) == 'list'){
    baf <- vcf.map$BAF
    geno <- vcf.map$GT
    rownames(geno) <- rownames(baf) <- baf$Probe_Set_ID
  } else if(class(vcf.map) == 'big.matrix'){
    baf <- vcf.map
    ProbeSetID <- rownames(vcf.map)
    is.bm <- TRUE
  } else if (is.data.frame(vcf.map) | is.matrix(vcf.map)){
    baf <- vcf.map
    baf$Probe_Set_ID <- rownames(vcf.map)
  }
  
  
  # New fast method
  var.dat.m <- do.call(rbind, var.dat)
  var.dat.m$idx <- rep(names(var.dat), sapply(var.dat, nrow))
  var.dat.m$var.snps <- with(var.dat.m, log((num.snps / var)^-1)) # Max SNPs, Low Variance
  var.dat.m$probeset <- unlist(sapply(var.dat, rownames))
  if(is.bm){
    var.dat.m$init <- var.dat.m$probeset %in% ProbeSetID
  } else {
    var.dat.m$init <- var.dat.m$probeset %in% baf$Probe_Set_ID
  }
  var.dat.m <- var.dat.m[which(var.dat.m$init),]
  max.var.ids <- sapply(split(var.dat.m, var.dat.m$idx), function(i) {
    i[which.max(i$var.snps),]$probeset # Min variance (var), max coverage (num.snps), returns probeset (probeset)
  })
  baf.var <- baf[as.character(max.var.ids),,drop=FALSE]
  
  if(is.bm){
    baf.var <- as.data.frame(baf.var / 100)
    baf.var$Probe_Set_ID <- rownames(baf.var)
  }
  
  if(class(vcf.map)=='list'){
    geno.var <- geno[baf.var$Probe_Set_ID,]
    return(list("BAF"=baf.var,
                "GT"=geno.var,
                "Var"=NULL)) #fill in with index of variance data
  } else{
    return(baf.var)
  }
  
}