#' segmentDrift
#' @description Segments a matrix (D) given a set of genomic coordinates
#' stored in fdat using either CBS or PCF
#' @param segmenter either CBS or PCF (Default: PCF)
#' @param fdat Genomic data-frame where first two columns are "chrom" and "pos"
#' @param D Matrix corresponding to fdat containing distances to plot, samples in 
#' columns, rows are the genomic pos
#'
#' @return CNA object
#' @export
segmentDrift <- function(segmenter='PCF', fdat, D, kmin=5){
  if(any(colnames(fdat)[1:2] != c('chrom', 'pos'))){
    warning(paste0("Column names: ", paste(colnames(fdat)[1:2], collapse=","), 
                   " are not 'chrom' and 'pos' and will be replaced"))
    colnames(fdat)[1:2] <- c("chrom", "pos")
  }
  
  CNAo <- switch(segmenter,
                 "PCF"={
                   require(dplyr)
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
                   require(DNAcopy)
                   CNAo <- with(fdat, #[names(ra.lm$residuals),],
                                CNA(genomdat=D,
                                    chrom=as.factor(chrom),
                                    maploc=pos,
                                    data.type="logratio",
                                    sampleid=colnames(D)))
                   smoothed.CNAo <- smooth.CNA(CNAo)
                   seg.CNAo <- segment(smoothed.CNAo,alpha = 0.001, eta=0.05, verbose=1, min.width=5)
                   seg.CNAo
                 })
  return(CNAo)
}


#' bafDrift
#' @definition Calcualtes the amount of genetic drift in a sample
#' 
#' @param sample.mat 
#' @param debug should be set to FALSE and only changed when debugging
#' @return
#' @export
bafDrift <- function(sample.mat, debug=FALSE, ...){
  require(DNAcopy)
  data(snp6.dat)
  ## Get pairwise distance between loci
  M <- CCLid:::.normBAF(sample.mat)
  D.l <- list()
  
  ## Order based on genomic position
  match.idx <- match(rownames(M), snp6.dat$SNP$Probe_Set_ID)
  if(any(is.na(match.idx))) {
    rm.idx <- which(is.na(match.idx))
    M <- M[-rm.idx,]
    match.idx <- match.idx[-rm.idx]
  }
  M <- M[order(match.idx),]
  g.loci <- snp6.dat$SNP[which(snp6.dat$SNP$Probe_Set_ID %in% rownames(M)),]
  
  ## calculate distance
  if(ncol(M) > 10) stop(paste0("Too many samples being compared for drift: n=", ncol(M)))
  while(ncol(M) > 1){
    D <- apply(M, 2, function(m){
      M[,1] - m
    })
    D <- scale(D, scale=FALSE)
    D.l <- append(D.l, list(D))
    M <- M[,-1,drop=F]
  }
  names(D.l) <- colnames(sample.mat)[-ncol(sample.mat)]
  
  ## Segment (CBS/PCF) the difference
  cna.drift <- lapply(D.l, function(D, ...){
    seg.CNAo <- segmentDrift(fdat = as.data.frame(g.loci), D=D[,-1], ...)
    seg.CNAo$output <- .addSegSd(seg.CNAo)
    
    seg.drift <- CCLid:::.estimateDrift(seg.CNAo, z.cutoff=1:3)
    seg.CNAo$output <- seg.drift$seg
    class(seg.CNAo) <- 'CCLid'
    if(debug) CCLid:::plot.CCLid(seg.CNAo)
    return(list("frac"=seg.drift$frac,
                "cna.obj"=seg.CNAo))
  })
  
  
  return(list("frac"=lapply(cna.drift, function(i) i$frac),
              "cna.obj"=lapply(cna.drift, function(i) i$cna.obj)))
    
}

#' addSegSd
#' @description Adds segment SD to CBS segment objects
#'
#' @param seg.obj an object returned from DNAcopy::segment()
#'
#' @return
.addSegSd <- function(seg.obj){
  adj.segs <- lapply(split(seg.obj$output, f=seg.obj$output$ID), function(seg){
    seg.dat <- as.data.frame(seg.obj$data)
    seg.dat$chrom <- as.character(seg.dat$chrom)
    
    ## Loop through each segment and find SD of the raw data
    gr.dat <- makeGRangesFromDataFrame(seg.dat, seqnames.field='chrom', 
                                       start.field = c('maploc', 'pos'),
                                       end.field = c('maploc', 'pos'), 
                                       keep.extra.columns = TRUE)
    gr.seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = TRUE)
    seqlevelsStyle(gr.seg) <- seqlevelsStyle(gr.dat) <- 'UCSC'
    
    ov.idx <- findOverlaps(gr.dat, gr.seg)
    s.idx <- grep(unique(gr.seg$ID), colnames(mcols(gr.dat)), fixed = TRUE)
    sd.per.seg <- sapply(split(ov.idx, subjectHits(ov.idx)), function(ov.i){
      dat <- mcols(gr.dat[queryHits(ov.i),])[, s.idx]
      lim <- quantile(dat, probs=c(0.05, 0.95), na.rm=TRUE) ##winsorization
      dat[dat < lim[1] ] <- lim[1]
      dat[dat > lim[2] ] <- lim[2]
      round(mad(dat, na.rm = TRUE),3)
    })
    seg$seg.sd <- sd.per.seg
    
    return(seg)
  })
  
  return(do.call(rbind, adj.segs))
}

#' .estimateDrift
#' @description Estimates genetic drift given an SD adjusted DNAcopy::segment()
#' object.  Estimates this based on a t-statistic
#' 
#' @param seg.obj DNAcopy::segment() object
#' @param z.cutoff Default set to 1:4;  Finds segments higher than those t-statistic away from 0
#'
#' @return
#'
#' @examples
#' .estimateDrift(seg.CNAo, z.cutoff=1:3)
.estimateDrift <- function(seg.obj, ...){
  seg.gr <- makeGRangesFromDataFrame(seg.obj$output, keep.extra.columns = TRUE)
  drift.dat <- lapply(split(seg.gr, seg.gr$ID), function(seg, z.cutoff=c(1:4)){
    ##Calculate Z-score of each seg.mean
    seg.sd <- mean(rep(seg$seg.sd, (width(seg) / 1000000)), na.rm=TRUE)
    seg.z <- round((seg$seg.mean) / (seg.sd / sqrt(seg$num.mark)), 3) ## t
    #seg.z <- round((seg$seg.mean / seg.sd), 3) ## z
    seg.z2 <- round((seg$seg.mean / seg$seg.sd), 3) ## z
    seg.z4 <- round((seg$seg.mean) / (seg$seg.sd / sqrt(seg$num.mark)), 3) ## t
    seg$seg.z <- seg.z 
    
    frac.cnv <- round(width(seg) / sum(width(seg)),3)
    ## Create filter criteria
    diff.idx <- sapply(setNames(z.cutoff, z.cutoff), function(z){which(seg.z > z | seg.z < -z)})
    diff.sum <- sapply(setNames(diff.idx, z.cutoff), function(idx) {sum(frac.cnv[idx])})
    seg$t <- floor(abs(seg$seg.z))
    
    return(list("seg"=seg, "sum"=diff.sum))
  }, ...)

  ## Reform the DNAcopy segment object into original populated data structure
  seg.out <- unlist(as(sapply(drift.dat, function(i) i$seg), "GRangesList"))
  names(seg.out) <- NULL
  seg.out <- as.data.frame(seg.out)[,-c(4,5)]
  colnames(seg.out)[c(1:3)] <- c("chrom", "loc.start", "loc.end")
  seg.out <- seg.out[,c(colnames(seg.obj$output), "seg.z", "t")]
  
  return(list("frac"=sapply(drift.dat, function(i){ i$sum }),
              "seg"=seg.out))
}

#' sigDiffBaf
#' @returns GRanges object of significant BAF drift (z > 3)
#' for a cna.obj returned by bafDrift$cna.obj
#' 
#' @param each.sample 
#' @param sig.es 
#'
#' @return
#' @export
sigDiffBaf <- function(each.sample, sig.es=NULL){
  es <- each.sample$output
  sig.idx <- which(es$t >= 3)
  if(length(sig.idx) > 0){
    es <- es[sig.idx,] 
    sig.es <- lapply(split(es, es$ID), makeGRangesFromDataFrame, keep.extra.columns=TRUE)
  }
  return(sig.es)
}

#' .getDrift
#' @description Returns the genomic fraction of drift for significant
#' regions (z>3)
#' 
#' @param i 
#' @param idx 
#'
#' @return
#' @export
#'
#' @examples
.getDrift <- function(i, idx=1){
  if(length(i$frac) >= idx){
    ## Select the "z > 3" row from all baf-drift estimates
    ## idx = 1: Should be VCF compare to all matching cell lines
    ## idx = 2: Should be SNP cell line compared to all other matching SNP cell line
    if(is.list(i$frac)){
      delta <- i$frac[[idx]][3,,drop=FALSE]
    } else {
      delta <- i$frac[3,,drop=FALSE]
    }
    rownames(delta) <- gsub("RNA_", "", names(i$frac)[1])
  } else {
    delta <- NULL
  }
  if(!is.null(delta)) colnames(delta) <- gsub("_.*", "", colnames(delta))
  return(as.data.frame(delta))
}