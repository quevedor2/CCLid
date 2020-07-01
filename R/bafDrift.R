#' bafDrift
#' @description Calcualtes the amount of genetic drift in a sample
#' 
#' @param sample.mat Input sample matrix
#' @param centering Centering the metric on median, none, or mean
#' @param norm.baf Boolean to normalize BAF (Default =TRUE)
#' @param hom.filt.val Homozygous filtering threshold (Default = 0.07)
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#' @param ... Extra param
#' @param debug should be set to FALSE and only changed when debugging
#' @importFrom utils data
#' @importFrom stats na.omit
#' @importFrom stats median
#' 
#' @export
bafDrift <- function(sample.mat, debug=FALSE, centering='none', 
                     norm.baf=TRUE, hom.filt.val=0.07, snp6.dat, ...){
  #data(snp6.dat)
  ## Get pairwise distance between loci
  M <- if(norm.baf) .normBAF(sample.mat) else sample.mat
  hom.filt.idx <- (rowSums(M) <= (hom.filt.val * ncol(M)))
  if(any(na.omit(hom.filt.idx))) M <- M[-which(hom.filt.idx),]
  #M <- M[-which(apply(M, 1, median, na.rm=TRUE) == 0),]
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
    # Mx <- data.frame("index"=c(1:nrow(M)), "val1"=M[,1], "val2"=M[,2])
    # loessMod1 <- loess(val1 ~ index, data=Mx, span=0.10) # 10% smoothing span
    # loessMod2 <- loess(val2 ~ index, data=Mx, span=0.10) # 10% smoothing span
    # par(mfrow=c(1,1))
    # plot(M[,1], col=scales::alpha("blue", 0.3), pch=16)
    # points(M[,2], col=scales::alpha("red", 0.3), pch=16)
    # lines(predict(loessMod1, Mx[,'index', drop=FALSE]), col="blue")
    # lines(predict(loessMod2, Mx[,'index', drop=FALSE]), col="red")
    D <- apply(M, 2, function(m){
      M[,1] - m
    })
    D <- switch(centering,
                "median"={
                  colmed <- apply(D, 2, median, na.rm=TRUE)
                  D - matrix(rep(colmed, nrow(D)), byrow=TRUE, nrow=nrow(D))
                },
                "mean"=scale(D, scale=FALSE),
                D)
    D.l <- append(D.l, list(D))
    M <- M[,-1,drop=F]
  }
  names(D.l) <- colnames(sample.mat)[-ncol(sample.mat)]
  
  ## Segment (CBS/PCF) the difference
  cna.drift <- lapply(D.l, function(D, ...){
    # seg.CNAo <- CCLid::segmentDrift(fdat = as.data.frame(g.loci), D=D[,-1],
    #                                 rm.homo=FALSE, segmenter=segmenter)
    # seg.CNAo$output <- .addSegSd(seg.CNAo, winsorize.data=TRUE)
    # seg.drift <- .estimateDrift(seg.CNAo, z.cutoff=NULL)
    seg.CNAo <- CCLid::segmentDrift(fdat = as.data.frame(g.loci), D=D[,-1], 
                                    rm.homo=FALSE, ...)
    seg.CNAo$output <- .addSegSd(seg.CNAo, ...)
    seg.drift <- .estimateDrift(seg.CNAo, ...) #z.cutoff=c(0:5) [I suggest NULL]
    
    seg.CNAo$output <- seg.drift$seg
    class(seg.CNAo) <- 'CCLid'
    
    # pdf("~/test4.pdf")
    # ccl.id <- ccl.id # 'OVCAR-5'
    # meta.cclid <- meta.df[grep(paste0("^", ccl.id, "$"), meta.df$ID),]
    # scp.path <- "scp quever@192.168.198.99:"
    # path.tmp <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
    # cat(paste0(scp.path, file.path(path.tmp, "CCLE", "eacon", meta.cclid$CCLE, "ASCAT", "L2R", "*png "), paste0("CCLE_", meta.cclid$ID, ".png\n")))
    # cat(paste0(scp.path, file.path(path.tmp, "GDSC", "eacon", gsub(".cel", "", meta.cclid$GDSC, ignore.case=TRUE), "ASCAT", "L2R", "*png "), paste0("GDSC_", meta.cclid$ID, ".png\n")))
    # if(debug) CCLid:::plot.CCLid(seg.CNAo, min.z = 1) #
    # dev.off()
    return(list("frac"=seg.drift$frac,
                "cna.obj"=seg.CNAo))
  })
  
  
  return(list("frac"=lapply(cna.drift, function(i) i$frac),
              "cna.obj"=lapply(cna.drift, function(i) i$cna.obj)))
  
}