###########################
#### Preliminary steps ####
###########################
## Run before executing any of the following 3 analyses
benchmarkCCLid <- function(bench){
  library(VariantAnnotation)
  library(CCLid)
  
  ## Set dataset colors
  dataset.cols <- setNames(RColorBrewer::brewer.pal(6, "Dark2"),
                           c("GDSC", "CCLE", "gCSI", "CGP", "Pfizer"))
  
  ## Load in Ref mat file
  PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/CCLid/CCLid"
  analysis <- 'baf'
  ref.mat <- downloadRefCCL("BAF", saveDir = PDIR)
  format.dat <- formatRefMat(name="BAF", ref.mat=ref.mat, 
                             analysis='baf', bin.size=5e5)
  # var.dat10mb <- var.dat2
  # var.dat5mb <- var.dat2
  # var.dat500k <- var.dat
  # var.dat1mb <- var.dat2
  # var.dat <- var.dat2
  ref.mat <- format.dat$mat
  var.dat <- format.dat$var
  
  ref.mat.bkup <- ref.mat
  new.ids <- assignGrpIDs(ref.mat.bkup, meta.df)
  colnames(ref.mat) <- new.ids
}


#######################
#### Genetic Drift ####
#######################
## Measures the amount of genetic drift between cell lines using BAF
geneticDrift <- function(){
  vcf.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs/GDSC'
  data.type <- 'rna'
  
  vcf.file <- switch(data.type,
                      "rna"=setNames(c('EGAF00000660866.snpOut.vcf.gz'), 
                                     c("HT-29"))) #HT-29

  ## Load in the VCF of HT-29
  sample.mat <- mapVcf2Affy(file.path(vcf.dir, vcf.file))$BAF[,c('Probe_Set_ID', 'BAF'), drop=FALSE]
  rownames(sample.mat) <- sample.mat$Probe_Set_ID
  sample.mat <- as.data.frame(sample.mat[,-1,drop=FALSE])

  vcf.map.var <- mapVariantFeat(sample.mat, var.dat)
  vcf.to.use <- as.matrix(vcf.map.var[,1, drop=FALSE])
  rownames(vcf.to.use) <- vcf.map.var$Probe_Set_ID
  
  ## Append sample to matrix
  ov.idx <- overlapPos(comp = vcf.to.use, ref=ref.mat, 
                       mapping = 'probeset')
  x.mat <- cbind(vcf.to.use[ov.idx$comp], 
                 ref.mat[ov.idx$ref,])
  
  ## Calculate drift of Cell line with RNAseq with external control
  cl.idx <- grep(names(vcf.file), colnames(x.mat))
  em2.idx <- grep('_EM-2$', colnames(x.mat))
  x.drift <- bafDrift(sample.mat=x.mat[,c(cl.idx, 1, em2.idx)], segmenter='CBS')
  
  pdf(file.path(vcf.dir, paste0("drift_", names(vcf.file), ".pdf")), height = 4, width=5)
  plot(x.drift$cna.obj[[1]], low.sig.alpha=0, sample.size=70)
  dev.off()
  print(file.path(vcf.dir, paste0("drift_", names(vcf.file), ".pdf")))

  
  
}

#################################
#### Cell Line Decomposition ####
#################################
## Measures the detection of cell line decomposition
combinedCellLine <- function(){
  vcf.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs/GDSC'
  data.type <- 'rna'
  
  vcf.files <- switch(data.type,
                     "wes"=setNames(c('DU-145.sample_id.vcf'), c("DU-145")), #A549.sample_id.vcf [not A549]
                     "rna"=setNames(c('EGAF00000661931.snpOut.vcf.gz', 'EGAF00000660849.snpOut.vcf.gz'), 
                                    c("EM-2", "COLO 205"))) #EM-2 CaR-1
  
  ## Load in two VCFs (EM-2 and CaR-1) from a given technology
  vcf.maps <- lapply(vcf.files, function(f) {
    mapVcf2Affy(file.path(vcf.dir, f))$BAF[,c('Probe_Set_ID', 'BAF'), drop=FALSE]
  })
  sample.mat <- as.data.frame(Reduce(function(x,y) merge(x,y, by='Probe_Set_ID'), vcf.maps))
  rownames(sample.mat) <- sample.mat$Probe_Set_ID
  sample.mat <- sample.mat[,-1,drop=FALSE]
  colnames(sample.mat) <- names(vcf.files)
  
  vcf.map.var <- mapVariantFeat(sample.mat, var.dat)
  vcf.to.use <- vcf.map.var[,c(1,2)]
  rownames(vcf.to.use) <- vcf.map.var$Probe_Set_ID

  q <- seq(0, 1, by=0.1)
  all.deconv <- lapply(q, function(p){
    print(paste0("Proportion: ", p))
    prop <- c(p, 1-p)
    sample.x <- as.matrix(combineSamples('BAF', vcf.to.use, prop))
    colnames(sample.x) <- "X"
    
    ## Calculate samples with highest probability
    ov.idx <- overlapPos(comp = sample.x, ref=ref.mat, 
                         mapping = 'probeset')
    x.mat <- cbind(sample.x[ov.idx$comp], 
                   ref.mat[ov.idx$ref,])

    ## Logit model based on euclidean distance between samples
    x.dist <- similarityMatrix(x.mat, 'euclidean')
    x <- x.dist[,1]
    head(sort(x))
    D.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=meta.df)
    balanced <- balanceGrps(D.vals)
    D.model <- trainLogit(balanced, predictors=c('baf'))
    
    x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
    pred <- assemblePredDat(x.vals, known.class=FALSE)
    pred <- mkPredictions(pred, D.model)
    
    p.cols <- c('baf.fit', 'z', 'q')
    x.pred <- split(pred, pred$Var2)[[1]]
    x.pred$CL <- gsub("^.*?_", "", x.pred$Var1)
    
    ## ground-truth
    truth.cl <- lapply(setNames(names(vcf.files), names(vcf.files)), function(i){
      i.idx <- grep(paste0("_", i, "$"), x.pred$Var1)
      i.prob <- x.pred[i.idx, p.cols]
      rownames(i.prob) <- x.pred[i.idx,]$Var1
      return(i.prob)
    })
    
    ## predicted-matches
    sig.idx <- which(x.pred$baf.fit < 0.5)
    x <- x.pred[order(x.pred$z),]
    pred.prob <- x.pred[sig.idx, c(p.cols, "CL")]
    rownames(pred.prob) <- x.pred[sig.idx,]$Var1
    pred.cl <- split(pred.prob, pred.prob$CL)

    ## Deconvolute the samples using highest probability samples:
    A <- sample.x[ov.idx$comp,,drop=FALSE]
    M.t <- sapply(truth.cl, function(p.cl){
        combineSamples("BAF", ref.mat[ov.idx$ref, rownames(p.cl), drop=FALSE], 
                     rep(1/nrow(p.cl), nrow(p.cl)))
    })
    M.p <- sapply(pred.cl, function(p.cl){
      combineSamples("BAF", ref.mat[ov.idx$ref, rownames(p.cl), drop=FALSE], 
                     rep(1/nrow(p.cl), nrow(p.cl)))
    })
    
    # Decomposition with complete data
    #M1.mse <- .checkMse(A, M)
    M0 <- NNLM::nnmf(A, k = 0, check.k = FALSE, init=list(W0 = M.t));
    M0.deconv <- as.matrix(M0$H[,1] / colSums(M0$H)) # [0.9, 0.09, 0.01]
    rownames(M0.deconv) <- colnames(M.t)
    
    M1 <- NNLM::nnmf(A, k = 0, check.k = FALSE, init=list(W0 = M.p));
    M1.deconv <- as.matrix(M1$H[,1] / colSums(M1$H)) # [0.9, 0.09, 0.01]
    rownames(M1.deconv) <- colnames(M.p)
    
    return(list("truth"=truth.cl,
                "pred"=pred.cl,
                "nmf.truth"=M0.deconv,
                "nmf.pred"=M1.deconv))
  })
  names(all.deconv) <- as.character(q)
  
  ## Plotting....
  out.ids <- gsub(" ", "", paste(names(vcf.files), collapse="_"))
  pdf(file.path(vcf.dir, paste0("deconvolute_", out.ids, ".pdf")), height = 4, width=5)
  nmf.vectors <- c('nmf.truth', 'nmf.pred')
  lapply(setNames(nmf.vectors, nmf.vectors), function(i){
    sample.names <- colnames(sample.mat)
    nmf.d <- plyr::rbind.fill(lapply(all.deconv, function(dat) as.data.frame(t(dat[[i]]))))
    
    .reduceP <- function(all.deconv, type='truth', val='baf.fit'){
      lapply(all.deconv, function(x) { as.data.frame(t(sapply(x[[type]], function(i) mean(i[,val]))))})
    }
    nmf.prob <- switch(i,
                       "nmf.truth"=plyr::rbind.fill(.reduceP(all.deconv, 'truth', 'baf.fit')),
                       "nmf.pred"=plyr::rbind.fill(.reduceP(all.deconv, 'pred', 'baf.fit')))
    rownames(nmf.d) <- rownames(nmf.prob) <- q
    
    sample.cols <- setNames(RColorBrewer::brewer.pal(ncol(nmf.pred), "Set1"),
                            sample.names)
    
    split.screen(matrix(c(0, 1, 0, 0.7,
                          0, 1, 0.7, 1), nrow=2, byrow = TRUE))
    screen(2); par(mar=c(0, 4.1, 4.1, 2.1), xpd=FALSE);
    plot(0, type='n', xlim=c(0,1), ylim=c(0,1), las=1, 
         xaxt='n', ylab="P", yaxt='n')
    axis(side=2, at=c(0,1), labels=c(0.0, 1.0), las=1)
    for(cl in sample.names){
      ## Pred
      box.idx <- grep(cl, sample.names)
      
      rect(xleft = q - if(box.idx == 1) 0.03 else 0, 
           ybottom = rep(0, nrow(nmf.prob)), 
           xright = q + if(box.idx == 1) 0 else 0.03, 
           ytop = 1-nmf.prob[,cl], 
           border = NA, col=sample.cols[cl])
    }

    screen(1); par(mar=c(5.1, 4.1, 0.5, 2.1))
    plot(0, type='n', xlim=c(0,1), ylim=c(0,1), las=1, xaxt='n',
         xlab="Cellular proportion", ylab="Predicted proportion")
    axis(side = 1, at=seq(0, 1, by=0.1), labels=q, cex.axis=0.7)
    axis(side = 1, at=seq(1, 0, by=-0.1), line=1, tick = FALSE, labels=q, cex.axis=0.7)
    par(xpd=TRUE) # this is usually the default
    axis(side = 1, at= -0.1, labels=sample.names[1], line=0, tick=FALSE, cex.axis=0.8)
    axis(side = 1, at= -0.1, labels=sample.names[2], line=1, tick=FALSE, cex.axis=0.8)
    par(xpd=FALSE) # this is usually the default
    
    ## Add TRUTH lines
    abline(coef = c(0,1), col="grey")
    abline(coef = c(1,-1), col="grey")
    
    ## Add deconvolution lines
    for(cl in sample.names){
      ## Pred
      lines(x=rownames(nmf.d), y=nmf.d[,cl], 
            lty=switch(i, 'nmf.truth'=1, 'nmf.pred'=2), 
            col=sample.cols[cl], lwd=2)
      points(x=rownames(nmf.d), y=nmf.d[,cl], col=sample.cols[cl], 
             pch=switch(i, 'nmf.truth'=16, 'nmf.pred'=15))
    }
    close.screen(all.screens=TRUE)
  })
  dev.off()
  #print(file.path(vcf.dir, paste0("deconvolute_", out.ids, ".pdf")))
}

##########################
#### snpsCellIdentity ####
##########################
## Measures the number of SNPs needed to call cellular identity
snpsCellIdentity <- function(){
  vcf.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs/GDSC'
  data.type <- 'rna'
  
  vcf.file <- switch(data.type,
                     "wes"=setNames(c('DU-145.sample_id.vcf'), c("DU-145")), #A549.sample_id.vcf [not A549]
                     "rna"=setNames(c('EGAF00000661931.snpOut.vcf.gz'), c("EM-2"))) #EM-2
  
  ## Load in two VCFs (A549 and DU-145) from a given technology
  vcf.map <- mapVcf2Affy(file.path(vcf.dir, vcf.file))
  
  r <- c(seq(1, 0.1, by=-0.1), seq(0.1, 0.01, by=-0.01), seq(0.01, 0.001, by=-0.001))
  if(any(duplicated(as.character(r)))) r <- r[-which(duplicated(as.character(r)))]
  num.snps <- setNames(floor(length(var.dat) * r), r)
  
  
  frac.score <- lapply(r, function(frac){
    idx <- sort(sample(1:length(var.dat), size=frac * length(var.dat), replace = FALSE))
    vcf.map.var <- mapVariantFeat(vcf.map, var.dat[idx])
    vaf.to.map <- vcf.map.var
    
    ## Overlap the two VCFs to form a single matrix to combine
    ov.idx <- overlapPos(comp = vaf.to.map$BAF,
                         ref=ref.mat, mapping = 'probeset')
    x.mat <- cbind(vaf.to.map$BAF$BAF[ov.idx$comp], 
                   ref.mat[ov.idx$ref,])
    
    x.dist <- similarityMatrix(x.mat, 'euclidean')
    D.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=meta.df)
    balanced <- balanceGrps(D.vals)
    
    models <- trainLogit(balanced, predictors=c('baf'))
    x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
    pred <- assemblePredDat(x.vals, known.class=FALSE)
    pred <- mkPredictions(pred, models)
    
    p.cols <- c("Var1", "baf.fit", "z", "q")
    x.pred <- split(pred, pred$Var2)[[1]]
    # x.pred[order(x.pred$z),]
    # #tail(x.pred[order(x.pred$baf, decreasing = TRUE),])
    # cl.idx <- grep(paste0("_", names(vcf.file), "$"), x.pred$Var1)
    # cl.prob <- setNames(1-x.pred[cl.idx,]$baf.fit, x.pred[cl.idx,]$Var1)
    
    return(x.pred[,p.cols])
  })
  names(frac.score) <- r
  
  fits <- lapply(setNames(c('baf.fit', 'z'), c('baf', 'z')), function(f){
    fit <- Reduce(function(x,y) merge(x,y,by="Var1"), lapply(frac.score, function(i) i[,c('Var1', f)]))
    colnames(fit) <- c('ID', r)
    
    m.df <- melt(fit)
    m.df$variable <- num.snps[as.character(m.df$variable)]
    
    coi <- rep(scales::alpha("black", 0.5), nrow(m.df))
    coi[grep(paste0("CCLE_", names(vcf.file), "$"), m.df$ID)] <- scales::alpha(dataset.cols['CCLE'], 0.5)
    coi[grep(paste0("GDSC_", names(vcf.file), "$"), m.df$ID)] <- scales::alpha(dataset.cols['GDSC'], 0.5)
    m.df$clid <- coi
    return(m.df)
  })
  fit.style <- 'z'
  m.df <- fits[[fit.style]]
  if(fit.style == 'baf') m.df$value <- 1-m.df$value
  
  
  ## Plotting....
  pdf(file.path(vcf.dir, paste0(fit.style, "_", names(vcf.file), ".pdf")), height = 4, width=5)
  {
    boxplot(value ~ variable, data = m.df, col="#0000ff22", las=1, cex.axis=0.8,
            xlab="Number of SNPs", main=paste0(toupper(data.type), ": ", names(vcf.file)),
            ylab=switch(fit.style, "baf"="Probability", "z"="Z-statistic"),
            ylim=switch(fit.style, "baf"=c(0,1), "z"=c(-10, 5)),
            outline=FALSE, border=FALSE, xaxt='n')
    axis(side=1, at=c(1:length(num.snps)), labels=rep('', length(num.snps)), lwd.ticks = 0.5)
    axis(side = 1, at=seq(1, length(num.snps), by=2), labels = rev(num.snps[c(FALSE,TRUE)]), tick = TRUE, line=0, cex.axis=0.7)
    axis(side = 1, at=seq(2, length(num.snps), by=2), labels = rev(num.snps[c(TRUE,FALSE)]), tick = FALSE, line=1, cex.axis=0.7)
    if(fit.style=='z') abline(h = -3, lty=2)
    
    spl <- split(m.df, m.df$variable)
    ext.m.df <- do.call(rbind, lapply(spl, function(i){
      switch(fit.style,
             "baf"=i[i$value > quantile(i$value, 0.95),],
             "z"=i[i$value < quantile(i$value, 0.05),])
    }))
    beeswarm(value ~ variable, data = ext.m.df, method = 'swarm', corral = 'gutter',
             cex=0.8, pch = 16, pwcol=clid, add=TRUE)
  }
  dev.off()
  
  print(file.path(vcf.dir, paste0(fit.style, "_", names(vcf.file), ".pdf")))
}
