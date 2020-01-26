#' .vcf2AffyGT
#' @description Maps Genotypes called from VCFs (0/0, 0/1, 1/0, 1/1) to the Affymetrix
#' standard of 0, 1, 2 genotypes
#' 
#' @param vcf.gt VariantAnnotation object where the $GT column contains 1/1, 0/1, etc...
#'
#' @return
#'
#' @examples
#' .vcf2AffyGT(vcf.affy.gr)
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
#'
#' @return
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
#' @param vcf.gt 
#' @param affy.gt 
#'
#' @return
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
#' @return
#'
#' @examples
#' .normBAF(x=seq(0, 1, by=0.1), lower=F)
.normBAF <- function(x, lower=T){
  #x <- seq(0, 1, by=0.1)
  if(lower){
    nBAF <- (0.5 - abs(x - 0.5))
  } else {
    nBAF <- (0.5 + abs(x - 0.5))
  }
  
  nBAF
}

#' mapVcf2Affy
#' @description Takes a VCF files and reads it into VariantAnnotation package to overlap with
#' GRanges probeset dataset of snp6.dat.  It willc onvert the 0/0, 0/1, etc.. genotypes to 0,1,2.
#' Finally, it will reutnr the BAFs and normalized BAFs for each of the overlapping probesets.
#' 
#' @param vcfFile Absolute path to VCF file
#'
#' @return List composed of two dataframes, BAF and GT
#' @export
#'
#' @examples
#' vcfFile <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/denis_id/mutect_GDSC/EGAR00001252191_13305_1/EGAR00001252191_13305_1.vcf"
#' mapVcf2Affy(vcfFile)
mapVcf2Affy <- function(vcfFile){
  require(VariantAnnotation)
  message(paste0("Reading in VCF file (", basename(vcfFile), "..."))
  vcf.gr <- readVcfAsVRanges(vcfFile)
  seqlevelsStyle(vcf.gr) <- 'UCSC'
  vcf.gr <- sort(vcf.gr)
  
  ## Overlap VCF file with the snp6 SNP probesets
  ov.idx <- findOverlaps(vcf.gr, CCLid::snp6.dat$SNP)
  vcf.affy.gr <- vcf.gr[queryHits(ov.idx)]
  mcols(vcf.affy.gr) <- cbind(mcols(vcf.affy.gr), 
                              mcols(CCLid::snp6.dat$SNP[subjectHits(ov.idx),]))
  
  ## Add the Affymetrix Genotype Scores (0, 1, 2 from 0/0, 1/0, 0/1, and 1/1)
  affy.genotype <- .vcf2AffyGT(vcf.affy.gr)  ## e.g., Convert 0/0 -> 0
  vcf.affy.gr <- .revComp(vcf.affy.gr) ## Reverse complement negative strand alleles
  vcf.affy.gr$affyGT <- .fixGT(vcf.affy.gr, affy.genotype) ## Fix genotypes for negative strand
  if(any(vcf.affy.gr$affyGT == -1)) vcf.affy.gr <- vcf.affy.gr[which(vcf.affy.gr$affyGT != -1),]
  
  ## Calculate BAF
  vcf.affy.gr$BAF <- round(altDepth(vcf.affy.gr) / (refDepth(vcf.affy.gr) + altDepth(vcf.affy.gr)),2)
  flip.idx <- .fixGT(vcf.affy.gr, affy.genotype, ret.idx=T)
  vcf.affy.gr[flip.idx,]$BAF <- (1-vcf.affy.gr[flip.idx,]$BAF)
  #vcf.affy.gr$BAF <- round(altDepth(vcf.affy.gr) / totalDepth(vcf.affy.gr),2)
  vcf.affy.gr$nBAF <- .normBAF(vcf.affy.gr$BAF)
  
  vcf.baf.df <- mcols(vcf.affy.gr)[,c('Probe_Set_ID', 'BAF', 'nBAF')]
  vcf.gt.df <- mcols(vcf.affy.gr)[,c('Probe_Set_ID', 'GT', 'affyGT')]
  
  return(list("BAF"=vcf.baf.df,
              "GT"=vcf.gt.df))
}

#' mapVariantFeat
#' @description Selects only the features with the highest variance
#' per bin
#' 
#' @param vcf.map An object returned by mapVcf2Aff()
#' @param var.dat A list returned by formatRefMat()$var
#' @param slow.method 
#'
#' @return
#' @export
mapVariantFeat <- function(vcf.map, var.dat, slow.method=FALSE){
  if(class(vcf.map) == 'list'){
    baf <- vcf.map$BAF
    geno <- vcf.map$GT
    rownames(geno) <- rownames(baf) <- baf$Probe_Set_ID
  } else if(is.data.frame(vcf.map) | is.matrix(vcf.map)){
    baf <- vcf.map
    baf$Probe_Set_ID <- rownames(vcf.map)
  }
  
  
  if(!slow.method){
    # New fast method
    var.dat.m <- reshape::melt(var.dat)
    var.dat.m$variable <- gsub("^[0-9]*\\.", "", names(unlist(var.dat)))
    var.dat.m$init <- var.dat.m$variable %in% baf$Probe_Set_ID
    var.dat.m <- var.dat.m[which(var.dat.m$init),]
    max.var.ids <- sapply(split(var.dat.m, var.dat.m$L1), function(i) {
      i[which.min(i$value),]$variable
    })
    baf.var <- baf[as.character(max.var.ids),,drop=FALSE]
  } else {
    # Old slow method
    max.var <- lapply(var.dat, function(v){
      #max.var <- head(sort(v[names(v) %in% baf$Probe_Set_ID], decreasing = TRUE),1)
      max.var <- head(sort(v[names(v) %in% baf$Probe_Set_ID]),1)
      return(baf[names(max.var),])
    }) 
    baf.var <- do.call(rbind, max.var)
  }

  if(class(vcf.map)=='list'){
    geno.var <- geno[baf.var$Probe_Set_ID,]
    return(list("BAF"=baf.var,
                "GT"=geno.var,
                "Var"=NULL)) #fill in with index of variance data
  } else if(is.data.frame(vcf.map) | is.matrix(vcf.map)){
    return(baf.var)
  }
 
}
