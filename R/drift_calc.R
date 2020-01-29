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
    #seg.z <- round((seg$seg.mean) / (seg.sd / sqrt(seg$num.mark)), 3) ## t
    seg.z <- round((seg$seg.mean / seg.sd), 3) ## z
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