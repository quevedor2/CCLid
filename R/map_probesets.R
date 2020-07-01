#' .vcf2AffyGT
#' @description Maps Genotypes called from VCFs (0/0, 0/1, 1/0, 1/1) to the Affymetrix
#' standard of 0, 1, 2 genotypes
#' 
#' @param vcf.gt VariantAnnotation object where the $GT column contains 1/1, 0/1, etc...
#'
.vcf2AffyGT <- function(vcf.gt){
  vcf.affy.map <- c('1/1'=2, 
                    '0/1'=1,
                    '1/0'=1,
                    '0/0'=0)
  if(class(vcf.gt) == 'VRanges'){
    affy.gt <- vcf.affy.map[vcf.gt$GT]
  } else if (all(vcf.gt %in% names(vcf.gt))) {
    affy.gt <- vcf.affy.map[vcf.gt]
  } else {
    stop("Could not determine type of input data, must be vector of 0/0 genotypes or VariantAnnotation object")
  }
  if(any(is.na(affy.gt))) affy.gt[is.na(affy.gt)] <- -1
  
  return(affy.gt)
}

#' .revComp
#' @description Reverse complements the Alleles
#' 
#' @param vcf.gt VariantAnnotation ob ject
#' @param ret.idx  Boolean, if True, does not reverse complement any alleles
.revComp <- function(vcf.gt, ret.idx=FALSE){
  complementAllele <- c(A="T", T="A", C="G", G="C")
  neg.strand <- which(vcf.gt$Strand == '-')
  
  if(ret.idx){
    neg.strand
  } else {
    vcf.gt[neg.strand,]$Allele_A <- complementAllele[vcf.gt[neg.strand,]$Allele_A]
    vcf.gt[neg.strand,]$Allele_B <- complementAllele[vcf.gt[neg.strand,]$Allele_B]
    vcf.gt
  }
}

#' .fixGT
#' @description Switches the genotype to adjust for the reverse complement switch (0/0 -> 1/1)
#' 
#' @param vcf.gt VCF genotype (0,1,2)
#' @param affy.gt Affy6 genotype
#' @param ret.idx Returns the index of matching instead of genotype fix
#' @importFrom VariantAnnotation alt
.fixGT <- function(vcf.gt, affy.gt, ret.idx=FALSE){
  complementGenotype <- c('0'='2', 
                          '1'='1', 
                          '2'='0')
  rev.idx <- which(vcf.gt$Allele_A == alt(vcf.gt))
  if(ret.idx){
    rev.idx
  } else {
    affy.gt[rev.idx] <- complementGenotype[as.character(affy.gt[rev.idx])]
    as.integer(affy.gt)
  }
}

#' .normBAF
#' @description rescales the BAF from a 0-1 scale to a 0-0.5
#' 
#' @param x BAF vector
#' @param lower if TRUE, 0-0.5 range, ELSE, 0.5-1 range
#' 
.normBAF <- function(x, lower=T){
  #x <- seq(0, 1, by=0.1)
  x[x>1] <- 1; x[x<0] <- 0
  if(lower){
    nBAF <- (0.5 - abs(x - 0.5))
  } else {
    nBAF <- (0.5 + abs(x - 0.5))
  }
  
  nBAF
}

#' .jsonToGr
#' @description converts a JSON object into a VariantAnnotation object
#'
#' @param from.file Boolean to say whether it comes from a file or passed in object
#' @param json Json VCF format - OUTDATED
#' @importFrom VariantAnnotation VRanges
#' @importFrom IRanges IRanges
#' 
.jsonToGr <- function(json, from.file=TRUE){
  if(from.file){
    ad <- sapply(json, function(i) i$sampleinfo[[1]]$AD)
    names <- as.character(sapply(json, function(i) i$sampleinfo[[1]]$NAME))
    gt <- sapply(json, function(i) i$sampleinfo[[1]]$GT)
  } else {
    ad <- sapply(json$sampleinfo, function(i) i$AD)
    names <- sapply(json$sampleinfo, function(i) i$NAME)
    gt <- sapply(json$sampleinfo, function(i) i$GT)
  }
  null.idx <- sapply(ad, is.null)
  if(any(null.idx)) ad[which(null.idx)] <- '0,0'
  ad <- unlist(ad)
  vr <- VRanges(seqnames=as.character(if(from.file) sapply(json, function(i) i$chr) else json$chr), 
                ranges = IRanges(start=as.integer(if(from.file) sapply(json, function(i) i$pos) else json$pos), 
                                 end=as.integer(if(from.file) sapply(json, function(i) i$pos) else json$pos)),
                ref=as.character(if(from.file) sapply(json, function(i) i$ref) else json$ref),
                alt=as.character(if(from.file) sapply(json, function(i) i$alt) else json$alt),
                totalDepth=sapply(strsplit(ad, ','), function(i) sum(as.integer(i))),
                refDepth=as.integer(sapply(strsplit(ad, ','), function(i) i[[1]])),
                altDepth=as.integer(sapply(strsplit(ad, ','), function(i) i[[2]])),
                sampleNames = names)
  vr$GT <- as.character(gt)
  return(vr)
}

#' .filt
#' @description Removes probesets under a certain depth, especially useful
#' for RNAsequencing vcfs
#'
#' @param vcf list from mapVcf2Affy
#' @param min.depth Default=5
#'
#' @return A filtered list from mapVcf2Affy
.filt <- function(vcf, min.depth=5){
  rm.idx <- which(vcf[[1]]$tdepth <= min.depth)
  vcf[[1]] <- vcf[[1]][-rm.idx,,drop=FALSE]
  vcf[[2]] <- vcf[[2]][-rm.idx,,drop=FALSE]
  return(vcf)
}

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
