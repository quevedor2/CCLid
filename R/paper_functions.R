############################
#### drift_it.R Support ###
############################
### Functions for CN - BAF drift overlap
#' get Seg SD
#' @param gr.D genomic ranges of segment distances
#' @param gr.Draw GenomicRanges of raw probeset distances
#' @param winsorize.data Winsorize data boolean (Default=FALSE)
#' @param winsor winsorization threshold (Default=0.95)
#' @param n.scale Scaling factor for denominator
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom stats quantile
#' @importFrom stats sd
getSegSD <- function(gr.D, gr.Draw, winsorize.data=FALSE, winsor=0.95, n.scale=1){
  idx <- grep(paste0("^X?", unique(gr.D$ID), "$"), colnames(mcols(gr.Draw)))
  ov <- findOverlaps(gr.D, gr.Draw)
  
  std.err <- sapply(split(gr.Draw[subjectHits(ov),], queryHits(ov)), function(raw){
    dat <- mcols(raw)[,idx]
    if(winsorize.data){
      lim <- quantile(dat, c(winsor, 1-winsor), na.rm=TRUE)
      dat[dat > max(lim)] <- max(lim)
      dat[dat < min(lim)] <- min(lim)
    }
    sigma=sd(dat, na.rm=TRUE)
    n= sum(!is.na(dat))
    n=if(n.scale==0) n else n/(n * n.scale)
    return(sigma / sqrt(n))
  })
  return(round(std.err, 3))
}

#' addSegDat
#' @param CNAo DNAcopy object
#' @param ... extra params
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges mcols
#' 
addSegDat <- function(CNAo, ...){
  gr.Draw <- makeGRangesFromDataFrame(CNAo$data, start.field = "pos", end.field="pos", keep.extra.columns = TRUE)
  gr.seg <- makeGRangesFromDataFrame(CNAo$output, keep.extra.columns = TRUE)
  ids <- colnames(mcols(gr.Draw))
  std.errs <- sapply(split(gr.seg, gr.seg$ID)[ids], getSegSD, gr.Draw=gr.Draw, ...)
  
  CNAo$output$seg.sd <- as.numeric(unlist(std.errs[ids]))
  return(CNAo)
}

#' getBafDrifts
#' @description Gets the drift GRanges and genomic fraction of drift
#' for any pair of cell lines from different datasets given
#'
#' @param cl.pairs Vector of indices for cell line pairs in x.mat
#' @param x.mat Input matrix containing probeset BAF data
#' @param ref.ds Reference dataset (e.g. CCLE)
#' @param alt.ds Comparing dataset (e.g. GDSC)
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#' @param ... Extra param
#'
#' @return Drift object
#' @export
getBafDrifts <- function(cl.pairs, x.mat, ref.ds=NULL, alt.ds=NULL, 
                         snp6.dat, ...){
  ref.idx <- grep(paste0(ref.ds, "_"), colnames(x.mat)[cl.pairs])
  alt.idx <- grep(paste0(alt.ds, "_"), colnames(x.mat)[cl.pairs])
  all.idx <- c(ref.idx, alt.idx)
  
  if(length(all.idx) == 2){
    # sample.mat = x.mat[,cl.pairs[all.idx]]; debug = TRUE;
    # centering = centering; norm.baf = TRUE; segmenter='PCF', hom.filt.val=0
    # bafDrift(sample.mat = x.mat[,cl.pairs[all.idx]], debug = FALSE, segmenter='PCF',
    #          centering = centering, norm.baf = TRUE, hom.filt.val=0)
    bdf <- bafDrift(sample.mat = x.mat[,cl.pairs[all.idx]], hom.filt.val=0,
                    snp6.dat=snp6.dat, ...)
    #CCLid:::plot.CCLid(bdf$cna.obj[[1]], min.z=4)
    drift.score <- list("sig.gr"=bdf$cna.obj, #CCLid::sigDiffBaf(bdf$cna.obj[[1]]),
                        "frac"=bdf$frac[[1]])
  } else {
    drift.score <- list("sig.gr"=NULL, "frac"=NULL)
  }
  
  return(drift.score)
}

#' getCnDrifts
#' @description Gets the drift GRanges and genomic fraction of drift
#' for any pair of cell lines from different datasets given L2R data
#'
#' @param ref.l2r Matrix of L2R data of genomic bins by samples for reference dataset
#' @param alt.l2r Matrix of L2R data of genomic bins by samples for comparison dataset
#' @param fdat segment data 
#' @param seg.id element ID of ref.l2r that specifies the segment L2Rs
#' @param raw.id element ID of ref.l2r that specifies the raw probeset L2Rs
#' @param verbose Verbose
#' @param ... Extra param
#' @param cell.ids All cell line IDs to compare drift between
#' @param meta.df Cell line metadata, accessible from CCLid::ccl_table
#' 
#' @importFrom utils data
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom matrixStats rowDiffs
#' @importFrom stats median
#' @return CN drift object
#' @export
getCNDrifts <- function(ref.l2r, alt.l2r,fdat, seg.id, raw.id, 
                        cell.ids, verbose=TRUE, meta.df, ...){
  #data(meta.df)
  ## Index matching cell line pairs for the CN PSets
  ref.bin.ids <- assignGrpIDs(ref.l2r[[seg.id]], meta.df)
  alt.bin.ids <- assignGrpIDs(alt.l2r[[seg.id]], meta.df)
  alt.ref.idx <- data.frame("id"=as.character(cell.ids),
                            "ref"=as.integer(sapply(paste0("_", cell.ids, "$"), 
                                                    grep, x=ref.bin.ids)),
                            "alt"=as.integer(sapply(paste0("_", cell.ids, "$"), 
                                                    grep, x=alt.bin.ids)))
  na.idx <- apply(alt.ref.idx, 1, function(i)  any(is.na(i)))
  if(any(na.idx)) alt.ref.idx <- alt.ref.idx[-which(na.idx),]
  alt.ref.idx$id <- as.character(alt.ref.idx$id)
  rownames(alt.ref.idx) <- alt.ref.idx$id
  
  ## Create a distance betweeen L2R matrix:
  # idx <- c(grep("^NCI-H522$", alt.ref.idx$id), grep("^NB-1$", alt.ref.idx$id)) #83, 8
  # [idx,,drop=FALSE]
  cn.drift <- apply(alt.ref.idx, 1, function(ar.i, centering='none', quantnorm=FALSE, max.med=0.05){
    ref.idx = as.integer(ar.i['ref'])
    alt.idx = as.integer(ar.i['alt'])
    
    seg.data <- cbind(ref.l2r[[seg.id]][,ref.idx,drop=FALSE], alt.l2r[[seg.id]][,alt.idx,drop=FALSE])
    raw.data <- cbind(ref.l2r[[raw.id]][,ref.idx,drop=FALSE], alt.l2r[[raw.id]][,alt.idx,drop=FALSE])
    
    if(quantnorm){
      seg.data <- normalize.quantiles(seg.data)
      raw.data <- normalize.quantiles(raw.data)
    }
    
    D = rowDiffs(seg.data)
    Draw = rowDiffs(raw.data)
    
    D.l <- switch(centering,
                "median"={
                  D.med <- apply(D, 2, median, na.rm=TRUE)
                  Draw.med <- apply(Draw, 2, median, na.rm=TRUE)
                  list("seg"= D - matrix(rep(D.med, nrow(D)), byrow=TRUE, nrow=nrow(D)), 
                       "raw"= Draw - matrix(rep(Draw.med, nrow(Draw)), byrow=TRUE, nrow=nrow(Draw)))
                },
                "mean"={
                  list("seg"=scale(D, center=TRUE, scale=FALSE), 
                       "raw"=scale(Draw, center=TRUE, scale=FALSE))
                },
                "extreme"={
                  D.med <- apply(D, 2, median, na.rm=TRUE)
                  Draw.med <- apply(Draw, 2, median, na.rm=TRUE)
                  if(abs(D.med) >= max.med & abs(Draw.med) >= max.med){
                    if (verbose) print(paste0("Extreme difference in ", colnames(D)))
                    list("seg"= D - matrix(rep(D.med, nrow(D)), byrow=TRUE, nrow=nrow(D)), 
                         "raw"= Draw - matrix(rep(Draw.med, nrow(Draw)), byrow=TRUE, nrow=nrow(Draw)))
                  } else {
                    list("seg"=D, "raw"=Draw)
                  }
                },
                list("seg"=D, "raw"=Draw))
    
    return(D.l)
  }, ...) #centering='extreme'
  D = do.call(cbind, lapply(cn.drift, function(i) i$seg))
  Draw = do.call(cbind, lapply(cn.drift, function(i) i$raw))
  colnames(D) <- colnames(Draw) <- alt.ref.idx$id
  # save(D, Draw, alt.ref.idx, fdat, file="~/D2.rda")
  rm(ref.l2r, alt.l2r); gc()
  
  ## Segment and find discordant regions
  # CNAo <- CCLid::segmentDrift(fdat = fdat, D=D, segmenter=segmenter)
  # CNAo$data <- cbind(CNAo$data[,1:2], Draw)
  # sd.CNAo <- addSegDat(ids=alt.ref.idx$id[idx], CNAo=CNAo, winsorize.data=TRUE, n.scale=0.01)
  CNAo <- CCLid::segmentDrift(fdat = fdat, D=D, ...)
  sd.CNAo <- CNAo
  sd.CNAo$data <- cbind(sd.CNAo$data[,1:2], Draw)
  
  rm(D, Draw); gc()
  sd.CNAo$output <- .addSegSd(sd.CNAo, winsorize.data=TRUE)
  seg.drift <- .estimateDrift(sd.CNAo, data.type='lrr')
  sd.CNAo$output <- seg.drift$seg
  class(sd.CNAo) <- 'CCLid'
  
  # pdf("~/temp.pdf")
  # CCLid:::plot.CCLid(sd.CNAo, min.z=1)
  # dev.off()
  cn.drifts <- list("frac"=seg.drift$frac,
                    "cna.obj"=sd.CNAo)
  return(cn.drifts)
}

#' driftOverlapMetric
#' @description Calculate the amount of overlap between the CN and BAF inferred drift
#' GRanges objects
#'
#' @param gr.baf GRangesList of cell lines an their BAF-drifted regions
#' @param gr.cn GRangesList of cell lines an their CN-drifted regions
#' @param ov.frac overlap fraction cutoff: (Default: seq(0, 1, by=0.01))
#' @param cell.ids Vector of cell line names to check drift of
#' @param baf.z z threshold for BAF difference
#' @param cn.z z threshold for L2R difference
#' @param cn.gtruth if TRUE, isolates BAF for only CN drifted regions
#'
#' @importFrom IRanges findOverlapPairs
#' @importFrom GenomicRanges pintersect
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges mcols<-
#' @importFrom GenomicRanges width
#' @importFrom stats setNames
#' 
#' @return List containing the following elements:
#' "model" = non-linear least-square model fitted to concordance and sensitvity
#' "saturation"=Matrix of total saturation and concordance cutoff
#' "sens"=Matrix of sensitivity value for each concordance cutoff
#' "dat"=Matrix of concordance cutoff by cell lines
#' @export
driftOverlapMetric <- function(gr.baf, gr.cn, cell.ids, ov.frac=seq(0, 1, by=0.01),
                               baf.z=4, cn.z=2, cn.gtruth=FALSE){
  ## calculate genomic overlap metric with different concordance-thresholds
  ov.dat <- sapply(cell.ids, function(cid){
    if(!is.null(gr.cn[[cid]])){
      if(!is.null(gr.baf[[cid]])){
        ov.baf.cn <- findOverlapPairs(gr.baf[[cid]], gr.cn[[cid]])
        baf.cn <- pintersect(ov.baf.cn)
        mcols(baf.cn) <- NULL
        
        # x <- sapply(0:9, function(baf.z){
        #   ov.baf.cn@first$t >= baf.z
        # })
        # z <- apply(x, 2, function(i) (ov.baf.cn@second$t >= cn.z) == i)
        # apply(z, 2, table)
        # apply(z, 2, function(i){
        #   sum(width(baf.cn[which(i),])) / sum(width(baf.cn))
        # })
        
        baf.cn$baf <- ov.baf.cn@first$t >= baf.z
        baf.cn$cn <- ov.baf.cn@second$t >= cn.z
        
        if(cn.gtruth){
          ## Isolate for only CN drifted regions
          baf.cn <- baf.cn[which(baf.cn$cn),,drop=FALSE]
        }
        m.idx <- baf.cn$baf == baf.cn$cn
        conc.drift <- sum(width(baf.cn[which(m.idx),])) / sum(width(baf.cn))
        setNames(c(conc.drift, conc.drift > ov.frac), c("conc", ov.frac))
        
      } else {
        setNames(rep(NA, length(ov.frac)+1), c("conc", ov.frac))
      }
    } else {
      setNames(rep(NA, length(ov.frac)+1), c("conc", ov.frac))
    }
  })
  rm.idx <- which(is.na(colSums(ov.dat)))
  if(length(rm.idx) > 0) ov.dat <- ov.dat[,-rm.idx]
  
  ## Organize sensitivity and assign a non-linear least-squares model
  ov.df <- data.frame("y"=rev(rowSums(ov.dat[-1,], na.rm=TRUE) / ncol(ov.dat)),
                      "x"=as.numeric(rownames(ov.dat[-1,])))
  
  fit <- tryCatch({
    m<-with(ov.df, nls(y~ SSasymp(x, Asym, R0, lrc)))
    
    ## Calculate saturation points
    Asym_coef <- summary(m)$coefficients[1] ## horizontal asymptote
    R0_coef <- summary(m)$coefficients[2] ## response when x == 0
    lrc_coef <- summary(m)$coefficients[3] ## rate constant
    sat.values <- seq(0, 1, by=0.01)
    conc.sat <- sapply(setNames(sat.values, sat.values), function(saturation){
      satX <- Asym_coef * saturation
      -(log((satX-Asym_coef)/(R0_coef-Asym_coef))/exp(lrc_coef)) # Conc. threshold to reach X% sensitvity saturation
    })
    conc.sat <- as.matrix(round(conc.sat, 3))
    
    list("model"=m,
         "saturation"=conc.sat)
  }, error=function(e){
    list("model"=NA,
         "saturation"=NA)
  })
  
  
  return(append(fit, list("sens"=ov.df,
                          "dat"=ov.dat)))
}

### Functions for RNA (BAF) - SNP (BAF) drift overlap
#' getVcfDrifts
#' @description Check the drift of a given VCF file and all matching
#' cell line names based on the meta-data matched IDs
#'
#' @param vcfFile path to vcf file
#' @param meta.df Meta file linking RNA files to cell names
#' @param ref.dat Refrence data containing Reference matrix and variance
#' @param min.depth Minimum depth to consider for SNPs
#' @param dataset Either 'CCLE', 'GDSC', ro 'GNE'
#' @param centering Centering of data method to be passed into bafDrift() function
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#'
#' @return A list containing the fraction of genome drifted,
#' as well the significantly drifted regions CNAo
#' @export
getVcfDrifts <- function(vcfFile, ref.dat, meta.df,  
                         min.depth=5, centering='extreme',
                         dataset='GDSC', snp6.dat){
  vcf <- basename(vcfFile)
  cat(basename(vcf), "...\n")
  ## Load in VCF data and leftjoin to existing ref.mat
  vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, 
                        ref.mat=ref.dat$ref, min.depth=min.depth,
                        snp6.dat=snp6.dat)
  rna.idx <- switch(dataset,
                    "GDSC"=grep(gsub(".snpOut.*", "", vcf), meta.df$EGAF),
                    "CCLE"=grep(gsub(".snpOut.*", "", vcf), meta.df$SRR),
                    "GNE"=grep(gsub(".snpOut.*", "", vcf), meta.df$gCSI_RNA))  
  colnames(vcf.mat)[1] <- paste0("RNA_", meta.df[rna.idx, 'ID'])
  if(is.na(meta.df[rna.idx, 'ID']) & dataset=='GNE'){
    colnames(vcf.mat)[1] <- meta.df[rna.idx,'gCSI_RNA']
  }
  
  ## Identify matching cell line data to RNAseq
  ## Calculate drift of Cell line with RNAseq with external control
  match.idx <- grep(paste0("(^|_)", gsub("NCI-", ".*", meta.df[rna.idx,]$ID), "$"), colnames(vcf.mat))
  if(length(match.idx) > 1){
    x.drift <- bafDrift(sample.mat=vcf.mat[,match.idx, drop=FALSE], hom.filt.val=0,
                        norm.baf=TRUE, segmenter='PCF', centering=centering)
    # head(x.drift$cna.obj[[1]]$output, 40)
    # CCLid:::plot.CCLid(x.drift$cna.obj[[1]], min.z=2)
    summ=x.drift
    # ## Isolate siginificant different regions
    # sig.diff.gr <- lapply(x.drift$cna.obj, sigDiffBaf)
    # frac.cnt <- x.drift$frac
    # 
    # summ <- list("frac"=frac.cnt,
    #              "sig"=sig.diff.gr)
  } else {
    summ=NULL
  }
  return(summ)
}

#' summarizeFracDrift
#' @description Summarizes the $frac from the baf.drifts and cn.drifts
#' given a cutoff for cn or baf
#'
#' @param cn.drifts CN drifts - paper function
#' @param cn.z L2R z threshold
#' @param baf.drifts  BAF drifts
#' @param baf.z BAF z threshold
#' @param include.id include the IDs in the analysis
#'
#' @return list of BAF and CN summary frac
#' @export
summarizeFracDrift <- function(cn.drifts, cn.z, baf.drifts, 
                               baf.z, include.id=FALSE){
  ## Reduce CN fraction of the genome drifted to a matrix
  cn.frac <- plyr::rbind.fill(lapply(cn.drifts$frac, function(i) as.data.frame(t(i))))
  cn.ids <- names(cn.drifts$frac)
  cn.frac <- cn.frac[,order(as.integer(colnames(cn.frac)))]
  cn.frac[is.na(cn.frac)] <- 0
  rownames(cn.frac) <- cn.ids

  ## Reduce BAF fraction of the genome drifted to a matrix
  baf.frac <- lapply(baf.drifts, function(i) i$frac)
  baf.frac <- baf.frac[-which(sapply(baf.frac, is.null))]
  baf.ids <- names(baf.frac)
  baf.frac <- plyr::rbind.fill(lapply(baf.frac, function(i) {
    df.i <- as.data.frame(t(i))
    colnames(df.i) <- gsub("^D.", "", colnames(df.i))
    df.i
  }))
  baf.frac <- baf.frac[,order(as.integer(colnames(baf.frac)))]
  baf.frac[is.na(baf.frac)] <- 0
  rownames(baf.frac) <- baf.ids
  
  ## Create rowSums based on the cn.z and baf.z cutoffs
  cn.cols <- as.integer(colnames(cn.frac)) >= cn.z
  baf.cols <- as.integer(colnames(baf.frac)) >= baf.z
  cn.summ <- data.frame("nodrift"=rowSums(cn.frac[,which(!cn.cols),drop=FALSE]),
                        "drift"=rowSums(cn.frac[,which(cn.cols),drop=FALSE]))
  baf.summ <- data.frame("nodrift"=rowSums(baf.frac[,which(!baf.cols),drop=FALSE]),
                         "drift"=rowSums(baf.frac[,which(baf.cols),drop=FALSE]))
  if(include.id){ 
    cn.summ$ID <- cn.ids
    baf.summ$ID <- baf.ids
  }
  return(list("cn"=cn.summ,
              "baf"=baf.summ))
}

#' plotFracDrift
#'
#' @param summ.frac summ frac list to plot, contains $baf and $cn 
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom graphics abline
#' @importFrom graphics par
#' @importFrom graphics axis
#' @importFrom graphics polygon
#' @importFrom stats density
#' 
#' @export
plotFracDrift <- function(summ.frac){
  cn.baf.frac <- merge(summ.frac$baf, summ.frac$cn, by="ID", all=TRUE)
  
  cn.baf.d <- with(cn.baf.frac, drift.x - drift.y)
  min.val <- max(head(sort(cn.baf.d), 1))
  max.val <- min(tail(sort(cn.baf.d), 1))
  cn.baf.frac$max <- cn.baf.d >= max.val | cn.baf.d <= min.val
  # low.diff.idx <- head(order(abs(cn.baf.d)), 300)
  # low.ord <- cn.baf.frac[low.diff.idx,]
  # low.ord <- low.ord[order(rowSums(low.ord[,c('drift.x', 'drift.y')])),]
  # cn.baf.frac[match(tail(low.ord, 1)$ID, cn.baf.frac$ID),]$max <- TRUE
  
  par(mar=c(5.1, 4.1, 6, 6), xpd=FALSE)
  with(cn.baf.frac, plot(drift.x, drift.y, pch=16, col=scales::alpha("black", 0.6),
                         axes=FALSE,
                         xlab="Fraction drift (BAF)", ylab="Fraction drift (L2R)", 
                         xlim=c(0,1), ylim=c(0,1)))
  abline(coef=c(0,1), col="grey", lty=2)
  axis(side = 1, at = seq(0, 1, by=0.2))
  axis(side = 2, at = seq(0, 1, by=0.2))
  
  par(xpd=NA)
  dx <- density(cn.baf.frac$drift.x, na.rm=TRUE)
  polygon(x = c(min(dx$x), dx$x, 1), 
          y=c(1, scales::rescale(dx$y, to=c(1,1.1)), 1), col="black")
  
  dy <- density(cn.baf.frac$drift.y, na.rm=TRUE)
  polygon(y = c(min(dy$x), dy$x, 1), 
          x = c(1, scales::rescale(dy$y, to=c(1,1.1)), 1), col="black")
  
  with(cn.baf.frac[which(cn.baf.frac$max),],
       points(drift.x, drift.y, pch=16, col="red"))
  with(cn.baf.frac[which(cn.baf.frac$max),],
       text(drift.x + 0.01, drift.y, labels=ID, adj=0, cex=0.6))
}

############################
#### match_it.R Support ####
############################
#' splitToMNM
#' @description Reduces the prediction dataframe to samples with a predicted
#' Match and a ground truth of non-match (based on anontations)
#' @param p prediction data frame
#'
#' @return A subset of p
splitToMNM <-function(p){
  p.nm <- split(p, p$g.truth)[['NM']]
  p.m.nm <- split(p.nm, p.nm$baf.p.fit)[['M']]
  p.m.nm[p.m.nm == 'character(0)'] <- NA
  return(p.m.nm)
}

#' genErrBp
#' @description Takes the output of splitToMNM() and reduces it to a dataframe
#' of Errors and PCLs
#'  
#' @param p.m.nm  output matrix of SplitToMNM
#'
#' @return Dataframe of errs and pcls
genErrBp <- function(p.m.nm){
  errs <- sapply(strsplit(p.m.nm$cellosaurus, split="/"), function(i) paste(unique(gsub("\\[PCL\\]", "", i)), collapse=""))
  pcls <- sapply(strsplit(p.m.nm$cellosaurus, split="/"), function(i) any(grepl("PCL", i)))
  err.pcl <- data.frame("err"=errs, "pcl"=pcls, stringsAsFactors = FALSE)
  if(any(err.pcl == '')) err.pcl[err.pcl ==''] <- 'X'
  err.pcl$err <- factor(err.pcl$err, levels=c("X", "OI", "SO", "SS"))
  err.pcl$pcl <- factor(err.pcl$pcl, levels=c(FALSE, TRUE))
  return(err.pcl)
}

#' checkAgainst
#' @description Uses the assembled "prediction" matrix with CVCL ids attributed ot each
#' cell line to check whether the two cell line pairs are known to match in the 
#' Cellosaurus database
#'
#' @param mat A matrix containing "cvclA" and "cvclB" columns for cellosaurus IDs of cell lines
#' @param melt.cells Melted cellosaurus dataframe, accessible from CCLid::ccl_table
#' 
#' @importFrom utils data
#' @return Character vector of OI (originating in), SS (synonymous), SI (sample from), 
#' and PCL (problematic)
checkAgainst <- function(mat, melt.cells){
  #data(melt.cells)
  .getAcr <- function(A, B, fp){
    B.mat <- B == fp
    row.A <- apply(A == fp, 1, any)
    row.B <- apply(B.mat, 1, any)
    row.AB <- which(row.A == row.B)
    
    acr <- which(apply(B.mat[row.AB,,drop=FALSE], 2, any))
    return(if(length(acr) ==0) 0 else acr)
  }
  map <- apply(mat, 1, function(i){
    fpA <- fullpull(i['cvclA'], melt.cells)
    fpB <- fullpull(i['cvclB'], melt.cells)
    if(!is.null(fpA) & !is.null(fpB)){
      A=.getAcr(i['cvclA'], i['cvclB'], fpA)
      B=.getAcr(i['cvclB'], i['cvclA'], fpB)
      
      if(A==1) names(A) <- 'SS'; if(B==1) names(B) <- 'SS' 
      contA <- any(fpA[which(fpA$CVCL %in% i['cvclA']),]$PCL)
      contB <- any(fpB[which(fpB$CVCL %in% i['cvclB']),]$PCL)
      return(paste0(names(A), "/", names(B), 
                    if(contA | contB) "[PCL]"))
    } else {
      return("/")
    }
  })
  
  return(map)
}


###########################
#### drug_it.R Support ####
###########################
#' loadInPSets
#'
#' @param drug.pset PSets path from pharmacoGX downloads
#'
#' @export
loadInPSets <- function(drug.pset){
  # drug.pset <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/PSets'
  psets <- list("CCLE"=readRDS(file.path(drug.pset, "CCLE.rds")),
                "GDSC"=readRDS(file.path(drug.pset, "GDSC2.rds")),
                "GNE"=readRDS(file.path(drug.pset, "gCSI2.rds")))
  return(psets)
}

#' getCinScore
#'
#' @param psets PSets from pharmacoGX
#' @param cin.metric sum or mean to compute CIN70 score
#' @param cin70 CIN70 gene set, accessible from CCLid::ccl_table
#' 
#' @importFrom utils data
#' @importFrom PharmacoGx molecularProfiles
#' 
#' @return a CIN list
#' @export
getCinScore <- function(psets, cin.metric='sum', cin70){
  #data(cin70)
  
  rna <- lapply(psets, function(pset, mDataType='rnaseq'){
    mdat <- molecularProfiles(pset, mDataType)
    colnames(mdat) <- pset@molecularProfiles[[mDataType]]$cellid
    cin.idx <- unlist(sapply(cin70$ENS, grep, x=rownames(mdat)))
    mdat[cin.idx,]
  })
  
  cin <- switch(cin.metric,
         "sum"=lapply(rna, colSums),
         "mean"=lapply(rna, colMeans))
  
  return(cin)
}

#' getGeneExpr
#'
#' @param psets PSets from pharmacoGX
#' @param in.key for mapping gene IDs, type of Gene (Default=SYMBOL)
#' @param gene.id Gene ID to map to Ensembl IDs
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom PharmacoGx molecularProfiles
#' 
#' @export
getGeneExpr <- function(psets, gene.id, in.key='SYMBOL'){
  
  if(in.key != 'ENSEMBL'){
    ensids <- mapIds(org.Hs.eg.db, keys = gene.id, keytype = in.key, column="ENSEMBL")
  } else {
    ensids <- gene.id
  }
  
  
  rna <- lapply(psets, function(pset, mDataType='rnaseq'){
    mdat <- molecularProfiles(pset, mDataType)
    colnames(mdat) <- pset@molecularProfiles[[mDataType]]$cellid
    cin.idx <- match(ensids, gsub("\\..*", "", rownames(mdat)))
    na.idx <- is.na(cin.idx)
    # cin.idx <- unlist(sapply(ensids, grep, x=rownames(mdat)))
    if(any(na.idx)) {
      cin.idx <- cin.idx[-which(na.idx)]
      gene.id <- gene.id[-which(na.idx)]
    }
    mdat <- mdat[cin.idx,,drop=FALSE]
    rownames(mdat) <- gene.id
    mdat
  })
  
  return(rna)
}

#' corWithDrug
#' @description IN DEVELOPMENT
#'
#' @param dat.d Data D containg $tCIN
#' @param title titleof plot
#' @param text.thresh  text.threshold of 0.50
#' @param genes List of genes
#' @param col.idx column index
#'
#' @importFrom stats cor
#' @importFrom stats cor.test
#' @importFrom stats p.adjust
#' @importFrom stats setNames
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
corWithDrug <- function(dat.d, col.idx, title='', text.thresh=0.5, genes=NULL){
  dat.d.abc <- do.call(rbind, apply(dat.d[,-col.idx], 2, function(i){
    # plot(i, cn.d$tCIN)
    # plot(i, cn.d$drift)
    dat.d$tCIN
    df <- data.frame("n"=table((is.na(dat.d$tCIN) == FALSE) + (is.na(i) == FALSE))[['2']],
                     "cinR"=cor(dat.d$tCIN, i, use="complete.obs"),
                     "cinR.p"=tryCatch({cor.test(dat.d$tCIN, i, use="complete.obs")$p.value}, error=function(e){NA}),
                     "driftR"=cor(dat.d$drift, i, use="complete.obs"),
                     "driftR.p"=tryCatch({cor.test(dat.d$drift, i, use="complete.obs")$p.value}, error=function(e){NA}))
    cbind(df, t(sapply(genes, function(g){cor(dat.d[,g], i, use='complete.obs')})))
  }))
  dat.d.abc$cin.q <- p.adjust(dat.d.abc$cinR.p, method="fdr")
  dat.d.abc$drift.q <- p.adjust(dat.d.abc$driftR.p, method="fdr")
  dat.d.abc$n.scale <- scales::rescale(dat.d.abc$n, to=c(0.2, 2))
  # head(dat.d.abc[order(dat.d.abc$drift.q),],10)
  # head(dat.d.abc[order(dat.d.abc$cin.q),],10)
  
  ## Visualization
  with(dat.d.abc, plot(cinR, driftR, col=scales::alpha("black", 1), 
                       xlim=c(-0.5,0.5), ylim=c(-0.5,0.5), cex=n.scale, main=title))
  abline(h=0, v=0, col="black")
  
  q.thresh <- seq(0.1, 0.7, by=0.1)
  q.col <- setNames(rev(brewer.pal(length(q.thresh),"PuRd")), q.thresh)
  used.idx <- c()
  for(i in q.thresh){
    idx <- which(dat.d.abc$cin.q  < i | dat.d.abc$drift.q < i)
    uidx <- setdiff(idx, used.idx)
    sig.dat.d <- dat.d.abc[uidx,]
    with(sig.dat.d, points(cinR, driftR, col=scales::alpha(q.col[as.character(i)], 0.7), pch=16, cex=n.scale))
    if(i < text.thresh & nrow(sig.dat.d) >= 1){
      with(sig.dat.d, text(cinR+0.02, driftR, labels=rownames(sig.dat.d), col=q.col[as.character(i)], adj=0, cex=0.8))
    }
    used.idx <- c(used.idx, uidx)
  }
  
  legend("topright", col=q.col, legend = paste0("q <= ", names(q.col)), pch=16, cex = 0.8)
  return(dat.d.abc)
}


########################
#### PLTK Functions ####
########################
#' cnTools: get Genes in TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' @param genome.build Genome build, only supports 'hg19'
#'
#' @description Gets the genes from UCSC hg19 TxDb knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom IRanges elementNROWS
#' @importFrom GenomicRanges granges
#' @return A Granges object containing strand-specific genes with EntrezIDs
getGenes <- function(genome.build="hg19"){
  switch(genome.build,
         hg19={ 
           package <- TxDb.Hsapiens.UCSC.hg19.knownGene 
          },
         stop("genome must be hg18, hg19, or hg38"))
  
  genes0 <- genes(package)
  idx <- rep(seq_along(genes0), elementNROWS(genes0$gene_id))
  genes <- granges(genes0)[idx]
  genes$gene_id = unlist(genes0$gene_id)
  genes
}

#' cnTools: annotate GRanges segments
#'
#' @param cn.data [Data.frame]: A dataframe that can be converted to GRanges object easily, or a granges object
#' @param genes [GRanges]: A GRanges object of genes with gene_ids housing annotation data. Easiest as the output from getGenes()
#' @param mart [object]: A biomart object if you want to annotate missed genes with ensembl
#' @param use.mart [boolean]: If no mart is given, it will load in a biomart object
#'
#' @param out.key Gene out.key, default is SYMBOL
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectLength
#' @importFrom S4Vectors splitAsList 
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' 
#' @return Annotated GRanges object with gene ids for the input GRanges
#' @export
annotateSegments <- function(cn.data, genes, out.key="SYMBOL", mart=NULL, use.mart=FALSE){
  if(use.mart & is.null(mart)){
    mart <- useMart("ENSEMBL_MART_ENSEMBL")
    mart <- useDataset("hsapiens_gene_ensembl", mart)
  }
  
  if(class(cn.data) == 'GRanges') {
    gr0 <- cn.data
  } else {
    gr0 <- makeGRangesFromDataFrame(cn.data,keep.extra.columns=TRUE)
  }
  seqlevelsStyle(gr0) <- 'UCSC'
  olaps <- findOverlaps(genes, gr0, type="within")
  idx <- factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
  gr0$gene_ids <- splitAsList(genes$gene_id[queryHits(olaps)], idx)
  gr0$gene_ids <- lapply(gr0$gene_ids, function(input.id) {
    if(length(input.id) > 0){ 
      tryCatch({
        ens <- mapIds(org.Hs.eg.db,
                      keys=input.id,
                      column=out.key,
                      keytype="ENTREZID",
                      multiVals="first")
        
        if(!is.null(mart)){
          ens[which(is.na(ens))] <- getBM(
            mart=mart,
            attributes="external_gene_name",
            filters="entrezgene",
            values=names(ens[which(is.na(ens))]),
            uniqueRows=TRUE)
        }
        ens
      }, error=function(e){NULL})
    } else { NA }
  })
  return(gr0)
}
