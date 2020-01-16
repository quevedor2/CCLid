#' bafDrift
#' @definition Calcualtes the amount of genetic drift in a sample
#' 
#' @param sample.mat 
#'
#' @return
#' @export
bafDrift <- function(sample.mat){
  require(DNAcopy)
  ## Get pairwise distance between loci
  M <- .normBAF(sample.mat)
  D.l <- list()
  
  while(ncol(M) > 1){
    D <- apply(M, 2, function(m){
      M[,1] - m
    })
    D.l <- append(D.l, list(D))
    M <- M[,-1,drop=F]
  }
  names(D.l) <- colnames(sample.mat)[-ncol(sample.mat)]
  
  ## CBS the difference
  g.loci <- snp6.dat$SNP[match(rownames(sample.mat), snp6.dat$SNP$Probe_Set_ID),]
  
  lapply(D.l, function(D){
    CNAo <- CNA(genomdat=D[,-1,drop=FALSE], 
                chrom=as.character(seqnames(g.loci)),
                maploc=start(g.loci), 
                data.type="logratio",sampleid=colnames(D)[-1])
    smoothed.CNAo <- smooth.CNA(CNAo)
    seg.CNAo <- segment(smoothed.CNAo, verbose=1, alpha=0.001, eta=0.01, min.width=5)
    lapply(split(seg.CNAo$output, f=seg.CNAo$output$ID), function(seg){
      seg.gr <- makeGRangesFromDataFrame(seg, keep.extra.columns = T)
      
      ##Calculate Z-score of each seg.mean
      seg.sd <- sd(rep(seg.gr$seg.mean, (width(seg.gr) / 1000000)))
      seg.mu <- mean(rep(seg.gr$seg.mean, (width(seg.gr) / 1000000)))
      seg.z <- (seg.gr$seg.mean - seg.mu) / seg.sd
      
      frac.cnv <- round(width(seg.gr) / sum(width(seg.gr)),3)
      ## Create filter criteria
      diff.idx <- which(seg.z > 2 | seg.z < -2)
      #small.seg.idx <- which(frac.cnv <= quantile(frac.cnv, 0.1))
      #data.frame("Fraction"=x)
      sum(frac.cnv[diff.idx])
    })
    plot(seg.CNAo, plot.type="w")
  })
    
}