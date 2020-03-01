############################
#### drift_it.R Support ####
############################
### Functions for CN - BAF drift overlap
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

addSegDat <- function(ids, CNAo, ...){
  gr.Draw <- makeGRangesFromDataFrame(CNAo$data, start.field = "pos", end.field="pos", keep.extra.columns = TRUE)
  gr.seg <- makeGRangesFromDataFrame(CNAo$output, keep.extra.columns = TRUE)
  ids <- colnames(mcols(gr.Draw))
  std.errs <- sapply(split(gr.seg, gr.seg$ID)[ids], CCLid:::getSegSD, gr.Draw=gr.Draw, ...)
  
  CNAo$output$seg.sd <- as.numeric(unlist(std.errs[ids]))
  return(CNAo)
}

#' getBafDrifts
#' @description Gets the drift GRanges and genomic fraction of drift
#' for any pair of cell lines from different datasets given
#' @param cl.pairs Vector of indices for cell line pairs in x.mat
#' @param x.mat Input matrix containing probeset BAF data
#' @param ref.ds Reference dataset (e.g. CCLE)
#' @param alt.ds Comparing dataset (e.g. GDSC)
#'
#' @return Drift object
#' @export
getBafDrifts <- function(cl.pairs, x.mat, ref.ds=NULL, alt.ds=NULL, ...){
  ref.idx <- grep(paste0(ref.ds, "_"), colnames(x.mat)[cl.pairs])
  alt.idx <- grep(paste0(alt.ds, "_"), colnames(x.mat)[cl.pairs])
  all.idx <- c(ref.idx, alt.idx)
  
  if(length(all.idx) == 2){
    # sample.mat = x.mat[,cl.pairs[all.idx]]; debug = TRUE;
    # centering = centering; norm.baf = TRUE; segmenter='PCF', hom.filt.val=0
    # bafDrift(sample.mat = x.mat[,cl.pairs[all.idx]], debug = FALSE, segmenter='PCF',
    #          centering = centering, norm.baf = TRUE, hom.filt.val=0)
    bdf <- bafDrift(sample.mat = x.mat[,cl.pairs[all.idx]], hom.filt.val=0, ...)
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
#' @param ref.l2r Matrix of L2R data of genomic bins by samples for reference dataset
#' @param alt.l2r Matrix of L2R data of genomic bins by samples for comparison dataset
#' @param cell.ids All cell line IDs to compare drift between
#' @param segmenter Segmentation algorithm, either 'PCF' (fast) or 'CBS' (slow)
#'
#' @return CN drift object
#' @export
getCNDrifts <- function(ref.l2r, alt.l2r,fdat, seg.id, raw.id, cell.ids, ...){
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
      require(preprocessCore)
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
                    print(paste0("Extreme difference in ", colnames(D)))
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
  rm(ref.l2r, alt.l2r, bins); gc()
  
  ## Segment and find discordant regions
  # CNAo <- CCLid::segmentDrift(fdat = fdat, D=D, segmenter=segmenter)
  # CNAo$data <- cbind(CNAo$data[,1:2], Draw)
  # sd.CNAo <- addSegDat(ids=alt.ref.idx$id[idx], CNAo=CNAo, winsorize.data=TRUE, n.scale=0.01)
  CNAo <- CCLid::segmentDrift(fdat = fdat, D=D, ...)
  sd.CNAo <- CNAo
  sd.CNAo$data <- cbind(sd.CNAo$data[,1:2], Draw)
  
  rm(D, Draw); gc()
  sd.CNAo$output <- CCLid:::.addSegSd(sd.CNAo, winsorize.data=TRUE)
  seg.drift <- CCLid:::.estimateDrift(sd.CNAo, data.type='lrr')
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
#'
#' @return List containing the following elements:
#' "model" = non-linear least-square model fitted to concordance and sensitvity
#' "saturation"=Matrix of total saturation and concordance cutoff
#' "sens"=Matrix of sensitivity value for each concordance cutoff
#' "dat"=Matrix of concordance cutoff by cell lines
#' @export
driftOverlapMetric <- function(gr.baf, gr.cn, cell.ids, ov.frac=seq(0, 1, by=0.01),
                               baf.z=4, cn.z=2, cn.gtruth=FALSE){
  require(GenomicRanges)
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
#' @param vcfFile path to vcf file
#' @param rna.meta.df Meta file linking RNA files to cell names
#'
#' @return A list containing the fraction of genome drifted,
#' as well the significantly drifted regions CNAo
#' @export
getVcfDrifts <- function(vcfFile, ref.dat, rna.meta.df, ...){
  vcf <- basename(vcfFile)
  cat(basename(vcf), "...\n")
  ## Load in VCF data and leftjoin to existing ref.mat
  vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, 
                        ref.mat=ref.dat$ref, ...)
  rna.idx <- switch(dataset,
                    "GDSC"=grep(gsub(".snpOut.*", "", vcf), rna.meta.df$EGAF),
                    "CCLE"=grep(gsub(".snpOut.*", "", vcf), rna.meta.df$SRR))  ## ADJUST THE GREP
  colnames(vcf.mat)[1] <- paste0("RNA_", rna.meta.df[rna.idx, 'ID'])
  
  ## Identify matching cell line data to RNAseq
  ## Calculate drift of Cell line with RNAseq with external control
  match.idx <- grep(paste0("_", gsub("NCI-", ".*", rna.meta.df[rna.idx,]$ID), "$"), colnames(vcf.mat))
  if(length(match.idx) > 1){
    x.drift <- bafDrift(sample.mat=vcf.mat[,match.idx, drop=FALSE], hom.filt.val=0,
                        norm.baf=TRUE, centering='median', segmenter='PCF')
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

#' readinRnaFileMapping
#' @description Map the RNA files to the SNP files
#' using hardcoded metadata
#' @return
#' @export
readinRnaFileMapping <- function(){
  require(CCLid)
  data(rna.meta.df)
  data(meta.df)
  
  all.meta.df <- merge(rna.meta.df, meta.df, by="ID", all.x=TRUE)
  return(all.meta.df)
}

#' summarizeFracDrift
#' @description Summarizes the $frac from the baf.drifts and cn.drifts
#' given a cutoff for cn or baf
#' @param cn.drifts 
#' @param cn.z 
#' @param baf.drifts 
#' @param baf.z 
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
#' @param summ.frac 
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
                         xlab="Fraction drift (BAF)", ylab="Fraction drift(L2R)", 
                         xlim=c(0,1), ylim=c(0,1)))
  abline(coef=c(0,1), col="grey", lty=2)
  axis(side = 1, at = seq(0, 1, by=0.2))
  axis(side = 2, at = seq(0, 1, by=0.2))
  
  par(xpd=NA)
  dx <- density(cn.baf.frac$drift.x)
  polygon(x = c(min(dx$x), dx$x, 1), 
          y=c(1, scales::rescale(dx$y, to=c(1,1.1)), 1), col="black")
  
  dy <- density(cn.baf.frac$drift.y)
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
#' @param mat A matrix containing "cvclA" and "cvclB" columns for cellosaurus IDs of cell lines
#'
#' @return Character vector of OI (originating in), SS (synonymous), SI (sample from), 
#' and PCL (problematic)
checkAgainst <- function(mat){
  require(Rcellosaurus)
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
#' @param drug.pset 
#' @export
loadInPSets <- function(drug.pset){
  # drug.pset <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/PSets'
  psets <- list("CCLE"=readRDS(file.path(drug.pset, "CCLE.rds")),
                "GDSC"=readRDS(file.path(drug.pset, "GDSC2.rds")),
                "GNE"=readRDS(file.path(drug.pset, "gCSI2.rds")))
  return(psets)
}

#' getCinScore
#' @param psets 
#' @param cin.metric 
#' @return a CIN list
#' @export
getCinScore <- function(psets, cin.metric='sum'){
  require(PharmacoGx)
  data(cin70)
  
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