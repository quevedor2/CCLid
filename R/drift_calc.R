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
  M <- CCLid:::.normBAF(sample.mat)
  D.l <- list()
  
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
  
  ## CBS the difference
  g.loci <- snp6.dat$SNP[match(rownames(sample.mat), snp6.dat$SNP$Probe_Set_ID),]
  
  cna.drift <- lapply(D.l, function(D){
    CNAo <- CNA(genomdat=D[,-1,drop=FALSE], 
                chrom=as.factor(seqnames(g.loci)),
                maploc=start(g.loci), 
                data.type="logratio",sampleid=colnames(D)[-1])
    smoothed.CNAo <- smooth.CNA(CNAo)
    seg.CNAo <- segment(smoothed.CNAo, verbose=1, alpha=0.01, eta=0.05, min.width=5)
    seg.CNAo$output <- .addSegSd(seg.CNAo)
    
    seg.drift <- .estimateDrift(seg.CNAo, z.cutoff=1:3)
    seg.CNAo$output <- seg.drift$seg
    class(seg.CNAo) <- 'CCLid'
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
                                       start.field = 'maploc', end.field = 'maploc', 
                                       keep.extra.columns = TRUE)
    gr.seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = TRUE)
    seqlevelsStyle(gr.seg) <- seqlevelsStyle(gr.dat) <- 'UCSC'
    
    ov.idx <- findOverlaps(gr.dat, gr.seg)
    s.idx <- grep(unique(gr.seg$ID), colnames(mcols(gr.dat)))
    sd.per.seg <- sapply(split(ov.idx, subjectHits(ov.idx)), function(ov.i){
      round(sd(mcols(gr.dat[queryHits(ov.i),])[, s.idx], na.rm = TRUE),3)
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
    seg.sd <- mean(rep(seg$seg.sd, (width(seg) / 1000000)))
    seg.z <- round((seg$seg.mean) / (seg.sd / sqrt(seg$num.mark)), 3)
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
  seg.out <- as.data.frame(seg.out)[,-c(4,5)]
  colnames(seg.out)[c(1:3)] <- c("chrom", "loc.start", "loc.end")
  seg.out <- seg.out[,c(colnames(seg.obj$output), "seg.z", "t")]
  
  return(list("frac"=sapply(drift.dat, function(i){ i$sum }),
              "seg"=seg.out))
}

#' plot.CCLid
#' plot() function for CCLid adjusted DNAcopy segment objects
#' 
#' @param obj 
#'
#' @return
#'
#' @examples
plot.CCLid <- function (obj, sample.size=50, low.sig.alpha=0.01, hi.sig.alpha=0.2) {
  require(scales)
  chroms <- paste0("chr", c(1:22, "X", "Y"))
  chr.data <- split(obj$data, obj$data$chrom)
  chr.seg <- split(obj$output, obj$output$chrom)
  samples <- colnames(obj$data)[-c(1:2)]
  chr.cols <- c('black', 'green')
  seg.col <- 'orange'
  sig.col <- 'red'
  
  uylim <- max(abs(obj$data[, -(1:2)]), na.rm = TRUE)
  ylim <- c(-uylim, uylim)

  ## Segment yb chromosome
  for(s in samples){
    par(mfrow=c(1,length(chroms)+2), mar=c(5.1, 0, 4.1, 0))
    plot(0, type='n', axes=FALSE, xlab='')
    plot(0, type='n', axes=FALSE, xlab='')
    
    for(chr in chroms){
      chr.col <- chr.cols[(match(chr, chroms) %% 2) + 1]
      sample.dat <- as.data.frame(chr.data[[chr]][,c('maploc', s), drop=FALSE])
      if(nrow(sample.dat) == 0) next
      
      size <- min(nrow(sample.dat), sample.size)
      sample.idx <- sample(1:nrow(sample.dat), size=size, replace=FALSE)
      sample.dat <- sample.dat[sample.idx,,drop=FALSE]
      
      plot(sample.dat, ylim=ylim, col=chr.col, pch=16, 
           xlim=c(1, max(sample.dat$maploc)), xlab='', 
           yaxt='n', xaxt='n', axes=FALSE, cex=0.6)
      abline(v = 1, lwd=0.5, col='grey')
      axis(side=1, at=median(sample.dat$maploc, na.rm=TRUE), 
           labels=gsub("chr", "", chr),
           tick=FALSE, line=(match(chr, chroms) %% 2))
      abline(h = 0, col="grey", lty=1, lwd=1)
      s.chr.seg <- chr.seg[[chr]]
      
      if(match(chr, chroms) == 1){
        par(xpd=TRUE)
        axis(side = 2, at=seq(-0.5, 0.5, by=0.5), labels=seq(-0.5, 0.5, by=0.5), las=1)
        par(xpd=FALSE)
      }
      
      ## Identifies significant different regions
      s.chr.seg <- s.chr.seg[which(s.chr.seg$ID %in% s),,drop=FALSE]
      if(any(s.chr.seg$t > 0)){
        sig.chr.seg <- s.chr.seg[which(s.chr.seg$t > 0),]
        sig.chr.seg$alpha <- 0
        sig.chr.seg$alpha[sig.chr.seg$t < 3] <- low.sig.alpha
        sig.chr.seg$alpha[sig.chr.seg$t >= 3] <- hi.sig.alpha
        with(sig.chr.seg, rect(xleft = loc.start, ybottom = ylim[1], 
                               xright = loc.end, ytop = ylim[2],
                               border=NA, col = scales::alpha(sig.col, alpha)))
      }
      
      ## Adds the SD interval
      with(s.chr.seg, rect(xleft = loc.start, ybottom = seg.mean - seg.sd, 
                           xright = loc.end, ytop = seg.mean + seg.sd,
                           border=NA, col = scales::alpha(seg.col, 0.4)))
      ## Adds the seg.mean line
      with(s.chr.seg, segments(x0 = loc.start, y0 = seg.mean, 
                               x1 = loc.end, y1 = seg.mean, 
                               lwd=3, col=seg.col))
    }
  }
}