.blankGr <- function(){
  makeGRangesFromDataFrame(data.frame("chr"='chrZ', "start"=1,  "end"=1))
}

.grepNA <- function(pattern, x){
  idx <- grep(pattern, x)
  if(length(idx) > 0) idx else NA
}

.getChrLength <- function(){
  require(BSgenome.Hsapiens.UCSC.hg19)
  chr.lengths = seqlengths(Hsapiens)[1:24]
  chr.len.gr <- makeGRangesFromDataFrame(data.frame("chrom"=names(chr.lengths),
                                                    "loc.start"=rep(1, length(chr.lengths)),
                                                    "loc.end"=chr.lengths))
  chr.len.gr$cum.end <- cumsum(as.numeric(end(chr.len.gr)))
  chr.len.gr$cum.start <- chr.len.gr$cum.end - (end(chr.len.gr) -1)
  chr.len.gr$cum.mid <- with(chr.len.gr, (cum.start + ((cum.end - cum.start)/2)))
  return(chr.len.gr)
}

#' multiDriftPlot
#'
#' @param seg 
#' @param chr.size.gr 
#' @param ref.ds 
#' @param alt.ds 
#'
#' @return
#' @export
multiDriftPlot <- function(seg, chr.size.gr=NULL, 
                           ref.ds=NULL, alt.ds=NULL){
  if(is.null(ref.ds)) stop("Requires input of ref.ds (GDSC or CCLE)")
  if(is.null(alt.ds)) stop("Requires input of alt.ds (GDSC or CCLE)")
  
  if(is.null(chr.size.gr)) chr.size.gr <- .getChrLength()
  chrs <- gsub("^chr", "", as.character(seqnames(chr.size.gr)))
  
  # seg <- seg.sig[-null.idx][[2]]
  grl <- as(unlist(seg), "GRangesList")
  grl.idx <- setNames(c(.grepNA(paste0(ref.ds, ".*", alt.ds), names(grl)),
                        .grepNA(paste0("RNA.*", alt.ds), names(grl)),
                        .grepNA(paste0("RNA.*", ref.ds), names(grl))),
                      c(paste0(ref.ds, "/", alt.ds),
                        paste0("RNA/", alt.ds),
                        paste0("RNA/", ref.ds)))
  
  ## Initialize the plotting space
  plot(0, type='n', xlim=c(1, max(chr.size.gr$cum.end)), ylim=c(0, 3),
       xlab="Chr", yaxt='n', ylab='', axes=FALSE)
  axis(side=1, at=chr.size.gr$cum.mid[c(TRUE,FALSE)], 
       labels = chrs[c(TRUE,FALSE)], tick=FALSE, cex.axis=0.8)
  axis(side=1, at=chr.size.gr$cum.mid[c(FALSE,TRUE)], line=0.5,
       labels = chrs[c(FALSE, TRUE)], tick=FALSE, cex.axis=0.8)
  axis(side=2, at=seq(0.5, 2.5, by=1), labels=names(grl.idx), cex.axis=0.8, 
       las=1, tick=FALSE, line=-1)
  
  lapply(seq_along(grl.idx), function(i){
    ## Draw empty chromosomes
    with(chr.size.gr, rect(xleft = cum.start, ybottom = rep(i - 0.9, length(chr.size.gr)), 
                           xright = cum.end, ytop = rep(i - 0.1, length(chr.size.gr)), 
                           col = "gray95", border="black", lwd=1))
    
    ## populate with sig.diff regions
    if(!is.na(grl.idx[i])){
      gr0 <- grl[[grl.idx[i]]]
      ov <- findOverlaps(gr0, chr.size.gr)
      gr0$cum.start <- start(gr0[queryHits(ov)]) + chr.size.gr[subjectHits(ov)]$cum.start
      gr0$cum.end <- end(gr0[queryHits(ov)]) + chr.size.gr[subjectHits(ov)]$cum.start
      
      with(gr0, rect(xleft = cum.start, ybottom = rep(i - 0.89, length(gr0)), 
                     xright = cum.end, ytop = rep(i - 0.11, length(gr0)), 
                     col = "#fb6a4a", border = NA))
    }
  })
}

driftOverlap <- function(seg, ref.ds=NULL, alt.ds=NULL){
  if(is.null(ref.ds)) stop("Requires input of ref.ds (GDSC or CCLE)")
  if(is.null(alt.ds)) stop("Requires input of alt.ds (GDSC or CCLE)")
  
  # seg <- seg.sig[-null.idx][[1]]
  if(is.null(unlist(seg))) {
    na.mat <- matrix(rep(NA, 3), ncol=1)
    return(list(na.mat, na.mat, na.mat))
  }
  grl <- as(unlist(seg), "GRangesList")
  grl.idx <- setNames(c(.grepNA(paste0(ref.ds, "_.*", alt.ds), names(grl)),
                        .grepNA(paste0("RNA_.*", alt.ds), names(grl)),
                        .grepNA(paste0("RNA_.*", ref.ds), names(grl))),
                      c(paste0(ref.ds, "/", alt.ds),
                        paste0("RNA/", alt.ds),
                        paste0("RNA/", ref.ds)))
  
  cs <- combn(x=1:3, m=2)
  drift.ov <- apply(cs, 2, function(i){
    # i <- unlist(cs[,2])
    gr1 <- if(is.na(grl.idx[i[1]]))  .blankGr() else grl[[grl.idx[i[1]]]]
    gr2 <- if(is.na(grl.idx[i[2]]))  .blankGr() else grl[[grl.idx[i[2]]]]
    grI <- intersect(gr1, gr2)
    
    wI <- sum(width(grI))
    wG1 <- sum(width(setdiff(gr1, grI)))
    wG2 <- sum(width(setdiff(gr2, grI)))
    wEmpty <- 0
    
    gr1.chr <- as.character(seqnames(gr1))
    gr2.chr <- as.character(seqnames(gr2))
    if(any(gr1.chr == 'chrZ') & any(gr2.chr == 'chrZ')){
      wEmpty <- 1
      wI <- 0
    }
    
    w.frac <- as.matrix(round(c(wG1, wI, wG2, wEmpty) / sum(c(wG1, wI, wG2, wEmpty)),2))
    rownames(w.frac) <- c(names(grl.idx)[i[1]], 'intersect', names(grl.idx)[i[2]], 'no_drift')
    return(as.data.frame(w.frac))
  })
  return(drift.ov)
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