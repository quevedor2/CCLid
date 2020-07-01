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

#' addSegSd
#' @description Adds segment SD to CBS segment objects
#'
#' @param seg.obj an object returned from DNAcopy::segment()
#' @param winsor Winsorization threshold (Default = 0.95)
#' @param verbose Default = FALSE
#' @param ... Extra para
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges mcols
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom stats quantile
#' @importFrom stats sd
.addSegSd <- function(seg.obj, winsor=0.95, verbose=FALSE, ...){
  adj.segs <- lapply(split(seg.obj$output, f=seg.obj$output$ID), function(seg){
    if (verbose) print(paste0(unique(seg$ID), "..."))
    seg.dat <- as.data.frame(seg.obj$data)
    seg.dat$chrom <- as.character(seg.dat$chrom)
    
    ## Loop through each segment and find SD of the raw data
    gr.dat <- makeGRangesFromDataFrame(seg.dat, seqnames.field='chrom', 
                                       start.field = c('maploc', 'pos'),
                                       end.field = c('maploc', 'pos'), 
                                       keep.extra.columns = TRUE)
    gr.seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = TRUE)
    seqlevelsStyle(gr.seg) <- 'UCSC'
    seqlevelsStyle(gr.dat) <- 'UCSC'
    
    ov.idx <- findOverlaps(gr.dat, gr.seg)
    s.idx <- grep(paste0("^", unique(gr.seg$ID), "$"), colnames(mcols(gr.dat)))
    sd.per.seg <- sapply(split(ov.idx, subjectHits(ov.idx)), function(ov.i, winsorize.data=FALSE){
      dat <- mcols(gr.dat[queryHits(ov.i),])[, s.idx]
      if(winsorize.data){
        # if (verbose) print("Winsorizing")
        lim <- quantile(dat, probs=c(winsor, 1-winsor), na.rm=TRUE) ##winsorization
        dat[dat < min(lim) ] <- min(lim)
        dat[dat > max(lim) ] <- max(lim)
      }
      
      return(round(sd(dat, na.rm = TRUE),3))
      # t.dat <- tryCatch({
      #   t.test(na.omit(dat))
      # }, error=function(e){
      #   data.frame("statistic"=NA, "p.value"=NA, "stderr"=NA)
      # })
      # setNames(round(c(t.dat$statistic, t.dat$p.value, t.dat$stderr),5),
      #          c("t", "p", "seg.sd"))
    }, ...)
    seg$seg.sd <- sd.per.seg
    # seg <- cbind(seg, abs(t(sd.per.seg)))
    # seg$'seg.sd' <- 0
    
    return(seg)
  })
  
  return(do.call(rbind, adj.segs))
}


#' .compSegs
#' @description compare the difference between BAF profiles between pairwise
#' samples using a t-test of the raw probeset data
#'
#' @param seg.obj A seg obj
#' @param verbose Verbose (Default = FALSE)
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges mcols<-
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom IRanges findOverlapPairs
#' @importFrom GenomicRanges pintersect
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom stats t.test
#'
#' @return A list of seg objects
.compSegs <- function(seg.obj, verbose=FALSE){
  spl.segs <- split(seg.obj$output, f=seg.obj$output$ID)
  seg.dat <- as.data.frame(as.matrix(seg.obj$data))
  gr.dat <- makeGRangesFromDataFrame(seg.dat, seqnames.field='chrom', 
                                     start.field = c('maploc', 'pos'),
                                     end.field = c('maploc', 'pos'), 
                                     keep.extra.columns = TRUE)
  dat.m <- as.matrix(mcols(gr.dat))
  storage.mode(dat.m) <- 'numeric'
  mcols(gr.dat) <- dat.m
  seqlevelsStyle(gr.dat) <- 'UCSC'
  
  
  diff.segs <- lapply(spl.segs, function(seg0){
    ## Take one segment
    gr.seg0 <- makeGRangesFromDataFrame(seg0, keep.extra.columns = TRUE)
    seg0.id <- unique(gr.seg0$ID)
    
    diff.seg1 <- lapply(spl.segs, function(seg1){
      ## Take a second segment to compare against
      gr.seg1 <- makeGRangesFromDataFrame(seg1, keep.extra.columns = TRUE)
      seg1.id <- unique(gr.seg1$ID)
      if(seg0.id == seg1.id) return(NULL)
      #if (verbose)(paste0(seg0.id, " - ", seg1.id))
      
      ## Find all overlap/intersects between segments
      ov.segs <- findOverlapPairs(gr.seg0, gr.seg1)
      seg <- pintersect(ov.segs)
      mcols(seg) <- NULL
      
      ## Find the raw SNP probes that populate those intersect segments
      dat.seg.ov <- findOverlaps(gr.dat, seg)
      s0.idx <- grep(seg0.id, colnames(mcols(gr.dat)), fixed = TRUE)
      s1.idx <- grep(seg1.id, colnames(mcols(gr.dat)), fixed = TRUE)
      
      ## Test the difference between the raw probe data
      test.diff <- sapply(split(dat.seg.ov, subjectHits(dat.seg.ov)), function(ov.i){
        dat <- mcols(gr.dat[queryHits(ov.i),])[, c(s0.idx, s1.idx)]
        dat.t <- t.test(dat[,1], dat[,2])
        return(round(c(dat.t$statistic, dat.t$p.value),3))
      })
      
      ## Append new metadata to the intersect/overlap between the segments
      new.meta <- cbind(mcols(ov.segs@first)[,c("ID", "seg.mean"),drop=FALSE],
                        mcols(ov.segs@second)[,c("ID", "seg.mean")],
                        t(test.diff))
      colnames(new.meta) <- c("refID", "seg.a", "ID", "seg.b", "seg.mean", "p")
      new.meta$t <- floor(abs(new.meta$seg.mean))
      mcols(seg) <- new.meta
      names(seg) <- NULL
      return(seg)
    })
    
    return(unlist(as(diff.seg1[-which(sapply(diff.seg1, is.null))], "GRangesList")))
  })

  return(diff.segs)
}


#' .estimateDrift
#' @description Estimates genetic drift given an SD adjusted DNAcopy::segment()
#' object.  Estimates this based on a t-statistic
#' 
#' @param ... Extra param
#' @param seg.obj DNAcopy::segment() object
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges width
#' @importFrom stats setNames
#' @importFrom methods as
#' 
.estimateDrift <- function(seg.obj, ...){
  seg.gr <- makeGRangesFromDataFrame(seg.obj$output, keep.extra.columns = TRUE)
  drift.dat <- lapply(split(seg.gr, seg.gr$ID), function(seg, z.cutoff=NULL, ...){
    ##Calculate Z-score of each seg.mean
    if('t' %in% colnames(mcols(seg))){
      seg$seg.z <- seg.z <- seg$t
    } else {
      # seg.sd <- mean(rep(seg$seg.sd, (width(seg) / 1000000)), na.rm=TRUE)
      seg.z <- round((seg$seg.mean / seg$seg.sd), 3) ## z
      seg$seg.z <- seg.z
    }
    frac.cnv <- round(width(seg) / sum(width(seg)),3)
    ## Create filter criteria
    if(is.null(z.cutoff)){
      seg <- .estimateZcutoff(seg, ...)
      diff.sum <- sapply(split(seg, seg$t), function(tseg) round(sum(width(tseg)) / sum(width(seg)),3))
    } else {
      diff.idx <- sapply(setNames(z.cutoff, z.cutoff), function(z){which(seg.z > z | seg.z < -z)})
      diff.sum <- sapply(setNames(diff.idx, z.cutoff), function(idx) {sum(frac.cnv[idx])})
      seg$seg.diff <- NA
      seg$t <- as.integer(cut(abs(seg$seg.z), breaks = z.cutoff, right = FALSE)) - 1
    }
    # seg$t <- floor(abs(seg$seg.z))
    
    return(list("seg"=seg, "sum"=diff.sum))
  }, ...)

  ## Reform the DNAcopy segment object into original populated data structure
  seg.out <- unlist(as(sapply(drift.dat, function(i) i$seg), "GRangesList"))
  names(seg.out) <- NULL
  seg.out <- as.data.frame(seg.out)[,-c(4,5)]
  colnames(seg.out)[c(1:3)] <- c("chrom", "loc.start", "loc.end")
  seg.out <- seg.out[,c(colnames(seg.obj$output), "seg.z", "t", "seg.diff")]
  
  return(list("frac"=sapply(drift.dat, function(i){ i$sum }),
              "seg"=seg.out))
}


#' .estimateZcutoff
#' @description Creates a theoreticla framework for difference
#' between BAFs, and then compares the observed z-differences
#' against the theoretical diff to look for differences
#'
#' @param seg an output dataframe from a CNA object
#' @param data.type Data type, supports only 'baf' at this time
#' @param verbose Verbose, default = FALSE
#'
#' @return Returned seg with $t and $seg.diff
.estimateZcutoff <- function(seg, data.type='baf', verbose=FALSE){
  ref.frac <- switch(data.type,
                     "baf"={
                       if (verbose) print("Estimating BAF diff significance from theoretical framework")
                       ref.frac <- sapply(c(1:5), function(tcn){
                         sapply(c(0:4), function(alt){
                           if(alt <= tcn) alt/tcn else 0
                         })
                       })
                       ref.frac <- round(unique(as.numeric(ref.frac)), 3)
                       ref.frac <- ref.frac[ref.frac <= 0.5]
                       (ref.frac)
                      },
                     "lrr"={
                       ref.frac <- seq(0, 1, by=0.1)
                       (ref.frac)
                     },
                     stop("data.type must be either baf or lrr"))

  ## Find all theoretical differences between BAFs in 100% purity state
  delta.frac <- abs(as.numeric(sapply(ref.frac, function(i) i -ref.frac)))
  delta.frac <- as.numeric(unique(as.character(round(sort(delta.frac),3))))
  
  ## Find which z-scores are greater than the theoretical SD-adjusted seg.mean difference
  theor.cutoff <- t(sapply(seg$seg.sd, function(s) round((delta.frac / s), 3)))
  t.mat <- (sweep(theor.cutoff, 1, as.matrix(abs(seg$seg.z))) <= 0)
  seg$t <-  rowSums(t.mat) - 1
  seg$seg.diff <- delta.frac[seg$t+1]
  # head(as.data.frame(seg), 30)
  return(seg)
}

#' sigDiffBaf
#' @returns GRanges object of significant BAF drift (z > 3)
#' for a cna.obj returned by bafDrift$cna.obj
#' 
#' @param each.sample Each sample
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' 
#' @export
sigDiffBaf <- function(each.sample){
  es <- each.sample$output
  sig.idx <- which(es$t >= 3)
  sig.es=NULL
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
#' @param i List object returned from bafDrift()
#' @param idx Index to return (1 = input compared to all matching cell lines, 2 = reference cell lines to all others)
#'
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
