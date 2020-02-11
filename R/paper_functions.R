############################
#### drift_it.R Support ####
############################
### Functions for CN - BAF drift overlap
getSegSD <- function(id, CNAo){
  idx <- grep(paste0("^X?", id, "$"), colnames(CNAo$data))
  D <- CNAo$data[,idx]
  
  number_of_chunks = ceiling(length(D) / 100)
  number_of_chunks=100
  sd.D <- sapply(split(seq_len(length(D)), 
                       cut(seq_len(length(D)), pretty(seq_len(length(D)), number_of_chunks))),
                 function(x) sd(D[x], na.rm=TRUE))
  sd.D <- mean(sd.D, na.rm=TRUE)
  return(sd.D)
}

addSegDat <- function(ids, CNAo){
  sds <- sapply(ids, getSegSD, CNAo=CNAo)
  CNAo$output$seg.sd <- rep(sds, table(CNAo$output$ID))
  seg.drift <- CCLid:::.estimateDrift(CNAo, z.cutoff=1:3)
  CNAo$output <- seg.drift$seg
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
    bdf <- bafDrift(x.mat[,cl.pairs[all.idx]], ...)
    #CCLid:::plot.CCLid(bdf$cna.obj[[1]])
    drift.score <- list("sig.gr"=bdf$cna.obj, #CCLid::sigDiffBaf(bdf$cna.obj[[1]]),
                        "frac"=bdf$frac[[1]][4,])
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
getCNDrifts <- function(ref.l2r, alt.l2r,fdat, cell.ids, segmenter='PCF'){
  ## Index matching cell line pairs for the CN PSets
  alt.bin.ids <- assignGrpIDs(ref.l2r, meta.df)
  ref.bin.ids <- assignGrpIDs(alt.l2r, meta.df)
  alt.ref.idx <- data.frame("id"=as.character(cell.ids),
                            "ref"=as.integer(sapply(paste0("_", cell.ids, "$"), 
                                                    grep, x=ref.bin.ids)),
                            "alt"=as.integer(sapply(paste0("_", cell.ids, "$"), 
                                                    grep, x=alt.bin.ids)))
  na.idx <- apply(alt.ref.idx, 1, function(i)  any(is.na(i)))
  if(any(na.idx)) alt.ref.idx <- alt.ref.idx[-which(na.idx),]
  alt.ref.idx$id <- as.character(alt.ref.idx$id)
  
  ## Create a distance betweeen L2R matrix:
  cn.drift <- apply(alt.ref.idx, 1, function(ar.i){
    ref.id = ref.bin.ids[ar.i['ref']]
    alt.id = alt.bin.ids[ar.i['alt']]
    
    ra.i <- as.data.frame(cbind(ref.l2r[,as.integer(ar.i['ref']), drop=FALSE],
                                alt.l2r[,as.integer(ar.i['alt']), drop=FALSE]))
    D = (ra.i[,1] - ra.i[,2])
    return(D)
  })
  rm(ref.l2r, alt.l2r); gc()
  
  ## Segment and find discordant regions
  CNAo <- switch(segmenter,
                 "PCF"={
                   CNdata <- with(fdat, cbind(as.factor(seg.seqnames), 
                                              seg.start, cn.drift))
                   colnames(CNdata) <- c("chrom", "pos", alt.ref.idx$id)
                   CNdata <- as.data.frame(CNdata)
                   pcf.dat <- copynumber::pcf(CNdata, pos.unit = "bp", kmin = 100, 
                                              gamma = 20, normalize = FALSE, 
                                              fast = TRUE, assembly = "hg19", 
                                              digits = 2, verbose = TRUE)
                   f <- factor(pcf.dat$sampleID, levels=alt.ref.idx$id)
                   pcf.dat <- pcf.dat[,c("sampleID", "chrom", "start.pos",
                                         "end.pos", "n.probes", "mean", "arm")]
                   colnames(pcf.dat) <- c("ID", "chrom", "loc.start", "loc.end", 
                                          "num.mark", "seg.mean", "arm")
                   
                   pcf.CNAo <- list("data"=CNdata,
                                    "output"=pcf.dat[order(f),],
                                    "segRows"=NULL,
                                    "call"=NULL)
                   pcf.CNAo
                 },
                 "CBS"={
                   CNAo <- with(fdat, #[names(ra.lm$residuals),],
                                CNA(genomdat=cn.drift,
                                    chrom=as.factor(seg.seqnames),
                                    maploc=seg.start,
                                    data.type="logratio",
                                    sampleid=alt.ref.idx$id))
                   smoothed.CNAo <- smooth.CNA(CNAo)
                   seg.CNAo <- segment(smoothed.CNAo,alpha = 0.001, eta=0.05, verbose=1, min.width=5)
                   seg.CNAo
                 })
  sd.CNAo <- CCLid:::addSegDat(ids=alt.ref.idx$id, CNAo=CNAo)
  seg.drift <- CCLid:::.estimateDrift(sd.CNAo, z.cutoff=1:4)
  sd.CNAo$output <- seg.drift$seg
  class(sd.CNAo) <- 'CCLid'
  # CCLid:::plot.CCLid(sd.CNAo)
  cn.drift <- list("frac"=seg.drift$frac,
                   "cna.obj"=sd.CNAo)
  return(cn.drift)
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
driftOverlapMetric <- function(gr.baf, gr.cn, cell.ids, ov.frac=seq(0, 1, by=0.01)){
  ## calculate genomic overlap metric with different concordance-thresholds
  ov.dat <- sapply(cell.ids, function(cid){
    if(! is.null(gr.cn[[cid]])){
      if(!is.null(gr.baf[[cid]])){
        ov.baf.cn <- findOverlapPairs(gr.baf[[cid]], gr.cn[[cid]])
        baf.cn <- pintersect(ov.baf.cn)
        mcols(baf.cn) <- NULL
        
        baf.cn$baf <- ov.baf.cn@first$t > 2
        baf.cn$cn <- ov.baf.cn@second$t > 2
        m.idx <- with(baf.cn, baf==cn)
        conc.drift <- sum(width(baf.cn[which(m.idx),])) / sum(width(baf.cn))
        
        setNames(conc.drift > ov.frac, ov.frac)
      } else {
        setNames(rep(0, length(ov.frac)), ov.frac)
      }
    } else {
      setNames(rep(NA, length(ov.frac)), ov.frac)
    }
    
    
  })
  
  ## Organize sensitivity and assign a non-linear least-squares model
  ov.df <- data.frame("y"=rev(rowSums(ov.dat, na.rm=TRUE) / ncol(ov.dat)),
                      "x"=as.numeric(rownames(ov.dat)))
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
  # rownames(conc.sat) <- rev(rownames(conc.sat))
  
  return(list("model"=m,
              "saturation"=conc.sat,
              "sens"=ov.df,
              "dat"=ov.dat))
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
getVcfDrifts <- function(vcfFile, ref.dat, rna.meta.df){
  vcf <- basename(vcfFile)
  cat(basename(vcf), "...\n")
  ## Load in VCF data and leftjoin to existing ref.mat
  vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, 
                        ref.mat=ref.dat$ref)
  rna.idx <- switch(dataset,
                    "GDSC"=grep(gsub(".snpOut.*", "", vcf), rna.meta.df$EGAF),
                    "CCLE"=grep(gsub(".snpOut.*", "", vcf), rna.meta.df$SRR))  ## ADJUST THE GREP
  colnames(vcf.mat)[1] <- paste0("RNA_", rna.meta.df[rna.idx, 'ID'])
  
  ## Identify matching cell line data to RNAseq
  ## Calculate drift of Cell line with RNAseq with external control
  match.idx <- grep(paste0("_", rna.meta.df[rna.idx,]$ID, "$"), colnames(vcf.mat))
  if(length(match.idx) > 0){
    x.drift <- bafDrift(sample.mat=vcf.mat[,match.idx], centering='mean', segmenter='PCF')
    
    ## Isolate siginificant different regions
    sig.diff.gr <- lapply(x.drift$cna.obj, sigDiffBaf)
    frac.cnt <- x.drift$frac
    
    summ <- list("frac"=frac.cnt,
                 "sig"=sig.diff.gr)
  } else {
    summ=NULL
  }
  return(summ)
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