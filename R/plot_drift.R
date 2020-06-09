#### Private Functions ####
###########################
#' .blankGr
.blankGr <- function(){
  makeGRangesFromDataFrame(data.frame("chr"='chrZ', "start"=1,  "end"=1))
}

#' .grepNA
#'
#' @param pattern pattern to grep
#' @param x vector to grep
.grepNA <- function(pattern, x){
  idx <- grep(pattern, x)
  if(length(idx) > 0) idx else NA
}

#' .getChrLength
.getChrLength <- function(){
  chr.lengths = seqlengths(Hsapiens)[1:24]
  chr.len.gr <- makeGRangesFromDataFrame(data.frame("chrom"=names(chr.lengths),
                                                    "loc.start"=rep(1, length(chr.lengths)),
                                                    "loc.end"=chr.lengths))
  chr.len.gr$cum.end <- cumsum(as.numeric(end(chr.len.gr)))
  chr.len.gr$cum.start <- chr.len.gr$cum.end - (end(chr.len.gr) -1)
  chr.len.gr$cum.mid <- (chr.len.gr$cum.start + ((chr.len.gr$cum.end - chr.len.gr$cum.start)/2))
  return(chr.len.gr)
}

#' .addCumPos
#'
#' @param dat Data containing $chrom, $pos, $loc.start, $loc.end
#' @param ref GRanges object containing $cum.start
#' @param dat.type Either "data" or "seg"
.addCumPos <- function(dat, ref, dat.type){
  m.row.idx <- match(as.character(dat$chrom), as.character(seqnames(ref)))
  if(dat.type=='data'){
    dat$cpos <- ref[m.row.idx,]$cum.start +  dat$pos - 1
    dat$chr.stat
  } else if(dat.type == 'seg'){
    dat$cloc.start <- ref[m.row.idx,]$cum.start +  dat$loc.start - 1
    dat$cloc.end <- ref[m.row.idx,]$cum.start +  dat$loc.end - 1
  }
  dat$chr.stat <- (m.row.idx %% 2) + 1
  return(dat)
}

#### Main Functions ####
########################
#' multiDriftPlot OUTDATED
#'
#' @param seg  Seg data
#' @param chr.size.gr NULL
#' @param ref.ds NULL
#' @param alt.ds NULL
multiDriftPlot <- function(seg, chr.size.gr=NULL, 
                           ref.ds=NULL, alt.ds=NULL){
  if(is.null(ref.ds)) stop("Requires input of ref.ds (GDSC or CCLE)")
  if(is.null(alt.ds)) stop("Requires input of alt.ds (GDSC or CCLE)")
  
  if(is.null(chr.size.gr)) chr.size.gr <- .getChrLength()
  chrs <- gsub("^chr", "", as.character(seqnames(chr.size.gr)))
  
  # seg <- seg.sig[-null.idx][[2]]
  grl <- as(unlist(seg), "GRangesList")
  grl.idx <- setNames(c(.grepNA(paste0(ref.ds, "_.*", alt.ds, "_"), names(grl)),
                        .grepNA(paste0("RNA_.*", alt.ds, "_"), names(grl)),
                        .grepNA(paste0("RNA_.*", ref.ds, "_"), names(grl))),
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
    rect(xleft = chr.size.gr$cum.start, ybottom = rep(i - 0.9, length(chr.size.gr)), 
         xright = chr.size.gr$cum.end, ytop = rep(i - 0.1, length(chr.size.gr)), 
         col = "gray95", border="black", lwd=1)
    
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

#' driftOverlap
#'
#' @param seg Seg data
#' @param ref.ds Reference dataset
#' @param alt.ds Alternate dataset
#'
#' @export
driftOverlap <- function(seg, ref.ds=NULL, alt.ds=NULL){
  if(is.null(ref.ds)) stop("Requires input of ref.ds (GDSC or CCLE)")
  if(is.null(alt.ds)) stop("Requires input of alt.ds (GDSC or CCLE)")
  
  # seg <- seg.sig[-null.idx][[1]]
  if(is.null(unlist(seg))) {
    na.mat <- matrix(rep(NA, 4), ncol=1)
    return(list(na.mat, na.mat, na.mat))
  }
  grl <- as(unlist(seg), "GRangesList")
  grl.idx <- setNames(c(.grepNA(paste0(ref.ds, "_.*", alt.ds), names(grl)),
                        .grepNA(paste0("RNA_.*", alt.ds), names(grl)),
                        .grepNA(paste0("RNA_.*", ref.ds), names(grl))),
                      c(paste0(ref.ds, "/", alt.ds),
                        paste0("RNA/", alt.ds),
                        paste0("RNA/", ref.ds)))
  
  drift.ov <- apply(combn(x=1:3, m=2), 2, function(i){
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
#' @param obj CCLid object returned from bafDrift()
#' @param sample.size Number of points to plot to allow for quicker plotting
#' @param low.sig.alpha low-confident drift regions alpha value for red
#' @param hi.sig.alpha high-confident drift regions alpha value for red
#' @param add.chr.sep Adds lines to delimit chromosomes
#' @param atype default=sd
#' @param add.points Adds points to the segments
#' @param min.z minimum z to report a drift change
#'
#' @examples
#' In the works...
plot.CCLid <- function (obj, sample.size=600, low.sig.alpha=0.01, 
                        hi.sig.alpha=0.2, add.chr.sep=TRUE, 
                        atype='sd', add.points=TRUE, min.z=3) {
  # low.sig.alpha=0.01
  # hi.sig.alpha=0.2
  # add.chr.sep=TRUE
  # sample.size=600
  # add.points=TRUE
  # atype='sd'
  chroms <- paste0("chr", c(1:22, "X", "Y"))
  if(any(grepl("(23)|(24)$", obj$data$chrom))){
    obj$data$chrom <- gsub("23$", "X", obj$data$chrom) %>%
      gsub("24$", "Y", .)
    obj$output$chrom <- gsub("23$", "X", obj$output$chrom) %>%
      gsub("24$", "Y", .)
  }
  if(any(!grepl("^chr", obj$data$chrom))){
    obj$data$chrom <- paste0("chr", obj$data$chrom)
    obj$output$chrom <- paste0("chr", obj$output$chrom)
  }
  
  chr.size.dat <- .getChrLength()
  
  chr.data <- .addCumPos(obj$data, chr.size.dat, dat.type='data')
  chr.seg <- .addCumPos(obj$output, chr.size.dat, dat.type='seg')
  samples <- colnames(obj$data)[-c(1:2)]
  chr.cols <- c('black', 'green')
  seg.col <- 'orange'
  sig.col <- 'red'
  
  ylim <- switch(atype,
                 "comp"=c(0,0.5),
                 c(-1,1))
  pos.col <- colnames(obj$data)[2]
  
  ## Segment yb chromosome
  for(s in samples){
    size <- min(nrow(chr.data), sample.size)
    sample.idx <- sample(1:nrow(chr.data), size=size, replace=FALSE)
    sample.dat <- chr.data[sample.idx, , drop=FALSE]
    sample.seg <- chr.seg[which(chr.seg$ID %in% s),]
    
    plot(x = sample.dat$cpos, y=sample.dat[,s], col=chr.cols[sample.dat$chr.stat],
         pch=16, ylim=ylim, xlim=c(1, max(chr.size.dat$cum.end)), xlab='', ylab=s,
         yaxt='n', xaxt='n', axes=FALSE, cex=0.6)
    abline(h = 0)
    if(add.chr.sep) abline(v = chr.size.dat$cum.end, lwd=0.5, col='grey')
    
    ## Adds chr labels
    axis(side=1, at = chr.size.dat$cum.mid[c(TRUE, FALSE)], 
         labels = gsub("^chr", "", as.character(seqnames(chr.size.dat)))[c(TRUE, FALSE)], 
         las=1,  tick=FALSE, lwd = 0, line=-1, cex.axis=0.6)
    axis(side=1, at = chr.size.dat$cum.mid[c(FALSE, TRUE)], 
         labels = gsub("^chr", "", as.character(seqnames(chr.size.dat)))[c(FALSE, TRUE)], 
         las=1,  tick=FALSE, lwd = 0, line=0, cex.axis=0.6)
    axis(side = 2, at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5), 
         las=1, cex.axis=0.8)
    
    ## Adds Signiifcant deviated regions
    if(any(na.omit(sample.seg$t) > 0)){
      sig.chr.seg <- sample.seg[which(sample.seg$t > 0),]
      sig.chr.seg$alpha <- 0
      sig.chr.seg$alpha[sig.chr.seg$t < min.z] <- low.sig.alpha
      sig.chr.seg$alpha[sig.chr.seg$t >= min.z] <- hi.sig.alpha
      with(sig.chr.seg, rect(xleft = cloc.start, ybottom = ylim[1], 
                             xright = cloc.end, ytop = ylim[2],
                             border=NA, col = scales::alpha(sig.col, alpha)))
    }
    
    ## Adds Segmean and SD intervals
    with(sample.seg, rect(xleft = cloc.start, ybottom = seg.mean - seg.sd, 
                       xright = cloc.end, ytop = seg.mean + seg.sd,
                       border=NA, col = scales::alpha(seg.col, 0.3)))
    with(sample.seg, segments(x0 = cloc.start, y0 = seg.mean, 
                           x1 = cloc.end, y1 = seg.mean, 
                           lwd=3, col=seg.col))
  }
}


#' plot.SimpleCCLid - OUTDATED
#'
#' @param obj CCLid object
#' @param add.chr.sep Adds chromosome separator
#' @param min.z Minimum Z for drift changes
plot.SimpleCCLid <- function(obj, add.chr.sep=TRUE, min.z=2){
  chroms <- paste0("chr", c(1:22, "X", "Y"))
  if(any(grepl("(23)|(24)$", obj$data$chrom))){
    obj$data$chrom <- gsub("23$", "X", obj$data$chrom) %>%
      gsub("24$", "Y", .)
    obj$output$chrom <- gsub("23$", "X", obj$output$chrom) %>%
      gsub("24$", "Y", .)
  }
  if(any(!grepl("^chr", obj$data$chrom))){
    obj$data$chrom <- paste0("chr", obj$data$chrom)
    obj$output$chrom <- paste0("chr", obj$output$chrom)
  }
  
  chr.size.dat <- .getChrLength()
  chr.seg <- .addCumPos(dat=obj$output, ref=chr.size.dat, dat.type='seg') #CCLid:::
  
  seg.col <- c('gray80', 'gray60')
  sig.col <- 'red'
  
  seg.spl <- split(chr.seg, chr.seg$ID)
  
  ### Viz begins!
  ## Create blank canvas to map the brilliance!
  plot(x = 0, type='n', xlab='', ylab='', 
       pch=16, ylim=c(0, length(seg.spl)), xlim=c(1, max(chr.size.dat$cum.end)), 
       yaxt='n', xaxt='n', axes=FALSE, cex=0.6)
  abline(h = 0)
  if(add.chr.sep) abline(v = chr.size.dat$cum.end, lwd=0.5, col='grey')
  
  ## Adds chr labels
  axis(side=1, at = chr.size.dat$cum.mid[c(TRUE, FALSE)], 
       labels = gsub("^chr", "", as.character(seqnames(chr.size.dat)))[c(TRUE, FALSE)], 
       las=1,  tick=FALSE, lwd = 0, line=-1, cex.axis=0.6)
  axis(side=1, at = chr.size.dat$cum.mid[c(FALSE, TRUE)], 
       labels = gsub("^chr", "", as.character(seqnames(chr.size.dat)))[c(FALSE, TRUE)], 
       las=1,  tick=FALSE, lwd = 0, line=0, cex.axis=0.6)
  axis(side = 2, at=(c(1:length(seg.spl)) - 0.5), labels=names(seg.spl), 
       las=1, cex.axis=0.8)
  
  for(s in names(seg.spl)){
    ssid <- seg.spl[[s]]
    ssid[is.na(ssid)] <- 0
    
    ssid$max.t <- ssid$t
    if(any(ssid$max.t > 5)) ssid[ssid$max.t > 5,]$max.t <- 5
    y.idx <- match(s, names(seg.spl))
    
    ## Add rectangles for the multiplots
    rect(xleft = ssid$cloc.start, ybottom = (y.idx - 0.8), 
         xright = ssid$cloc.end, ytop = (y.idx - 0.5 + ssid$max.t/10), 
         col=seg.col[ssid$chr.stat], border = NA)
    if(any(ssid$t >= min.z)){
      ss.sig <- ssid[which(ssid$t >= min.z),]
      rect(xleft = ss.sig$cloc.start, ybottom = (y.idx - 0.8), 
           xright = ss.sig$cloc.end, ytop = (y.idx - 0.5 + ss.sig$max.t/10), 
           col=sig.col, border = NA)
    }
  }
  
}


#' plot.multiObj OUTDATED
#'
#' @param OBJ A list of CCLid Objects
#' @param add.chr.sep Adds chromosome separator
#' @param min.z Z threshold to report drift change
#'
#' @export
plot.multiObj <- function(OBJ, add.chr.sep=TRUE, min.z=2){
  if(class(OBJ) == 'list'){
    par(mfrow=c(length(OBJ), 1), mar=c(3, 8, 3, 2))
    
    if(length(min.z) > 1){
      if(length(min.z) != length(OBJ)) stop("min.z length does not equal input length of CCLid objects")
      lapply(seq_along(OBJ), function(idx){
        plot.SimpleCCLid(OBJ[[idx]], min.z[idx], add.chr.sep=add.chr.sep)
      })
    } else {
      lapply(OBJ, plot.SimpleCCLid, min.z=min.z, add.chr.sep=add.chr.sep)
    }
  } else if(class(OBJ)=='CCLid'){
    par(mar=c(3, 8, 3, 2))
    plot.SimpleCCLid(OBJ, min.z=min.z, add.chr.sep=add.chr.sep)
  }
}
