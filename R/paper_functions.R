############################
#### drift_it.R Support ####
############################
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
getBafDrifts <- function(cl.pairs, x.mat, ref.ds=NULL, alt.ds=NULL){
  ref.idx <- grep(paste0(ref.ds, "_"), colnames(x.mat)[cl.pairs])
  alt.idx <- grep(paste0(alt.ds, "_"), colnames(x.mat)[cl.pairs])
  all.idx <- c(ref.idx, alt.idx)
  
  if(length(all.idx) == 2){
    bdf <- bafDrift(x.mat[,cl.pairs[all.idx]])
    #CCLid:::plot.CCLid(bdf$cna.obj[[1]])
    drift.score <- list("sig.gr"=CCLid::sigDiffBaf(bdf$cna.obj[[1]]),
                        "frac"=bdf$frac[[1]][3,])
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
getCNDrifts <- function(ref.l2r, alt.l2r, cell.ids, segmenter='PCF'){
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
                   CNdata <- with(featureData(bins[[alt.ds]])@data,
                                  cbind(as.factor(seg.seqnames), seg.start, cn.drift))
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
                   CNAo <- with(featureData(bins[[alt.ds]])@data, #[names(ra.lm$residuals),],
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