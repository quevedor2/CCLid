#' mapVcf2Affy
#' @description Takes a VCF files and reads it into VariantAnnotation package to overlap with
#' GRanges probeset dataset of snp6.dat.  It willc onvert the 0/0, 0/1, etc.. genotypes to 0,1,2.
#' Finally, it will reutnr the BAFs and normalized BAFs for each of the overlapping probesets.
#' 
#' @param vcfFile Absolute path to VCF file
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#' 
#' @importFrom rjson fromJSON
#' @importFrom VariantAnnotation readVcfAsVRanges
#' @importFrom VariantAnnotation refDepth
#' @importFrom VariantAnnotation altDepth
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges mcols<-
#' @importFrom S4Vectors subjectHits
#' 
#' @return List composed of two dataframes, BAF and GT
#' @export
#'
mapVcf2Affy <- function(vcfFile, snp6.dat){
  if(class(vcfFile) == 'list'){
    message(paste0("JSON data passed in..."))
    vcf.gr <- .jsonToGr(vcfFile, from.file=FALSE)
  } else if(grepl("\\.json$", vcfFile)){
    message(paste0("Reading in JSON file (", basename(vcfFile), "..."))
    json_data <- fromJSON(file=vcfFile)
    vcf.gr <- .jsonToGr(json_data, from.file=TRUE)
  } else if(grepl('\\.vcf', vcfFile)){
    message(paste0("Reading in VCF file (", basename(vcfFile), "..."))
    vcf.gr <- readVcfAsVRanges(vcfFile) # Memory: 700Mb -> 959Mb
  }
  seqlevelsStyle(vcf.gr) <- 'UCSC'
  vcf.gr <- sort(vcf.gr)
  
  ## Overlap VCF file with the snp6 SNP probesets
  ov.idx <- findOverlaps(vcf.gr, snp6.dat$SNP) # Memory: 1.6Gb
  vcf.affy.gr <- vcf.gr[queryHits(ov.idx)]
  mcols(vcf.affy.gr) <- cbind(mcols(vcf.affy.gr), 
                              mcols(snp6.dat$SNP[subjectHits(ov.idx),]))
  vcf.affy.gr$tdepth <- (refDepth(vcf.affy.gr) + altDepth(vcf.affy.gr))
  gc() # Memory: 1.4Gb
  
  ## Add the Affymetrix Genotype Scores (0, 1, 2 from 0/0, 1/0, 0/1, and 1/1)
  affy.genotype <- .vcf2AffyGT(vcf.affy.gr)  ## e.g., Convert 0/0 -> 0
  vcf.affy.gr <- .revComp(vcf.affy.gr) ## Reverse complement negative strand alleles
  vcf.affy.gr$affyGT <- .fixGT(vcf.affy.gr, affy.genotype) ## Fix genotypes for negative strand
  if(any(vcf.affy.gr$affyGT == -1)) vcf.affy.gr <- vcf.affy.gr[which(vcf.affy.gr$affyGT != -1),]
  
  ## Calculate BAF
  vcf.affy.gr$BAF <- round(altDepth(vcf.affy.gr) / (vcf.affy.gr$tdepth),2)
  flip.idx <- .fixGT(vcf.affy.gr, affy.genotype, ret.idx=T)
  vcf.affy.gr[flip.idx,]$BAF <- (1-vcf.affy.gr[flip.idx,]$BAF)
  #vcf.affy.gr$BAF <- round(altDepth(vcf.affy.gr) / totalDepth(vcf.affy.gr),2)
  vcf.affy.gr$nBAF <- .normBAF(vcf.affy.gr$BAF)
  gc() # Memory: 1.4Gb
  
  vcf.baf.df <- mcols(vcf.affy.gr)[,c('Probe_Set_ID', 'BAF', 'nBAF', 'tdepth')]
  vcf.gt.df <- mcols(vcf.affy.gr)[,c('Probe_Set_ID', 'GT', 'affyGT', 'tdepth')]
  
  return(list("BAF"=vcf.baf.df,
              "GT"=vcf.gt.df))
}