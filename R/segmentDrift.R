#' segmentDrift
#' @description Segments a matrix (D) given a set of genomic coordinates
#' stored in fdat using either CBS or PCF
#'
#' @param segmenter either CBS or PCF (Default: PCF)
#' @param fdat Genomic data-frame where first two columns are "chrom" and "pos"
#' @param D Matrix corresponding to fdat containing distances to plot, samples in 
#' columns, rows are the genomic pos
#' @param kmin Minimum number of SNPs to consider (default = 5)
#' @param rm.homo  Remove homozygous SNPs (Default = FALSE)
#' @importFrom stats median
#' @importFrom dplyr %>%
#' @import magrittr
#' @importFrom DNAcopy smooth.CNA
#' @importFrom DNAcopy segment
#' 
#' @return CNA object
#' @export
segmentDrift <- function(segmenter='PCF', fdat, D, kmin=5, rm.homo=FALSE){
  if(any(colnames(fdat)[1:2] != c('chrom', 'pos'))){
    warning(paste0("Column names: ", paste(colnames(fdat)[1:2], collapse=","), 
                   " are not 'chrom' and 'pos' and will be replaced"))
    colnames(fdat)[1:2] <- c("chrom", "pos")
  }
  
  if(rm.homo){
    med.val <- apply(D, 1, median, na.rm=TRUE)
    rm.idx <- which(med.val < 0.02) 
    fdat <- fdat[-rm.idx,,drop=FALSE]
    D <- D[-rm.idx,,drop=FALSE]
  }
  CNAo <- switch(segmenter,
                 "PCF"={
                   CNdata <- with(fdat, cbind(as.factor(chrom), pos, D))
                   CNdata <- as.data.frame(CNdata)
                   pcf.dat <- copynumber::pcf(CNdata, pos.unit = "bp", kmin = kmin, 
                                              gamma = 20, normalize = FALSE, 
                                              fast = TRUE, assembly = "hg19", 
                                              digits = 2, verbose = FALSE)
                   f <- factor(pcf.dat$sampleID, levels=colnames(D))
                   pcf.dat <- pcf.dat[,c("sampleID", "chrom", "start.pos",
                                         "end.pos", "n.probes", "mean", "arm")]
                   colnames(pcf.dat) <- c("ID", "chrom", "loc.start", "loc.end", 
                                          "num.mark", "seg.mean", "arm")
                   
                   colnames(CNdata)[1:2] <- c('chrom', 'pos')
                   CNdata$chrom <- gsub("^23$", "X", CNdata$chrom) %>% gsub("^24$", "Y", .) %>% gsub("^", "chr", .)
                   pcf.dat$chrom <- gsub("^23$", "X", pcf.dat$chrom) %>% gsub("^24$", "Y", .) %>% gsub("^", "chr", .)
                   
                   pcf.CNAo <- list("data"=CNdata,
                                    "output"=pcf.dat[order(f),],
                                    "segRows"=NULL,
                                    "call"=NULL)
                   pcf.CNAo
                 },
                 "CBS"={
                   CNAo <- with(fdat, #[names(ra.lm$residuals),],
                                DNAcopy::CNA(genomdat=D,
                                             chrom=as.factor(chrom),
                                             maploc=pos,
                                             data.type="logratio",
                                             sampleid=colnames(D)))
                   smoothed.CNAo <- DNAcopy::smooth.CNA(CNAo)
                   seg.CNAo <- DNAcopy::segment(smoothed.CNAo,alpha = 0.01, eta=0.05, verbose=1, min.width=5)
                   seg.CNAo
                 })
  return(CNAo)
}
