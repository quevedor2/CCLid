bencharkCCLid <- function(bench){
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
}


#################################
#### Cell Line Decomposition ####
#################################
## Measures the detection of cell line decomposition
combinedCellLine <- function(){
  ref.mat.bkup <- ref.mat
  new.ids <- assignGrpIDs(ref.mat.bkup, meta.df)
  colnames(ref.mat) <- new.ids
  
  vcf.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs/GDSC'
  data.type <- 'rna'
  
  vcf.files <- switch(data.type,
                     "wes"=setNames(c('DU-145.sample_id.vcf'), c("DU-145")), #A549.sample_id.vcf [not A549]
                     "rna"=setNames(c('EGAF00000661931.snpOut.vcf.gz', 'EGAF00000660849.snpOut.vcf.gz'), 
                                    c("EM-2", "COLO-205"))) #EM-2 CaR-1
  
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

  q <- seq(0, 1, by=0.05)
  all.deconv <- lapply(q, function(p){
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
    sig.idx <- which(x.pred$z < -3)
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
  
  lapply(all.deconv, function(i) i$pred)
  lapply(all.deconv, function(i) i$nmf.pred)
  lapply(all.deconv, function(i) i$nmf.truth)
  
  all.deconv <- as.data.frame(t(do.call(cbind, all.deconv)))
  all.deconv$prop <- q
    

  prob.mat <- as.data.frame(do.call(rbind, all.probs))
  prob.mat$prop_A <- q
  prob.mat$prop_B <- 1-q
  
  dir.create("~/CCLid")
  save(prob.mat, all.probs, file="~/CCLid/proportion_prob.rda")
  
  pdf("~/proportion.pdf")
  plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xaxt='n', xlab="Proportion", 
       ylab="P", las=1, main="Proportion of in silico mixed cell lines")
  cols <- setNames(c("#d73027", "#fc8d59", "#91bfdb", "#4575b4"),
                   colnames(prob.mat)[1:4])
  axis(side = 1, at=prob.mat$prop_A, labels=prob.mat$prop_A, line=0, tick = TRUE)
  axis(side = 1, at=prob.mat$prop_A, labels=prob.mat$prop_B, line=1, tick = FALSE)
  for(i in colnames(prob.mat[1:4])){
    points(x = prob.mat$prop_A, y=prob.mat[,i], pch=16, col=cols[i])
    lines(x = prob.mat$prop_A, y=prob.mat[,i], col=cols[i])
  }
  
  lines(x=all.deconv$prop, all.deconv$V1, lty=2, col="grey")
  lines(x=all.deconv$prop, all.deconv$V2, lty=2, col="grey")
  legend(x=0, y=0.9, fill=cols, legend=names(cols), box.lwd = 0)
  dev.off()
}

##########################
#### snpsCellIdentity ####
##########################
## Measures the number of SNPs needed to call cellular identity
snpsCellIdentity <- function(){
  ref.mat.bkup <- ref.mat
  new.ids <- assignGrpIDs(ref.mat.bkup, meta.df)
  colnames(ref.mat) <- new.ids
  
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
  fit.style <- 'baf'
  m.df <- fits[[fit.style]]
  if(fit.style == 'baf') m.df$value <- 1-m.df$value
  
  
  pdf(file.path(vcf.dir, paste0(fit.style, "_", names(vcf.file), ".pdf")), height = 4, width=6)
  boxplot(value ~ variable, data = m.df, col="#0000ff22", las=1, cex.axis=0.8,
          xlab="Number of SNPs", main=paste0(toupper(data.type), ": ", names(vcf.file)),
          ylab=switch(fit.style, "baf"="Probability", "z"="Z-statistic"),
          ylim=switch(fit.style, "baf"=c(0,1), "z"=c(-10, 5)),
          outline=FALSE, border=FALSE, xaxt='n')
  axis(side=1, at=c(1:length(num.snps)), labels=rep('', length(num.snps)), lwd.ticks = 0.5)
  axis(side = 1, at=seq(1, length(num.snps), by=2), labels = rev(num.snps[c(FALSE,TRUE)]), tick = TRUE, line=0, cex.axis=0.5)
  axis(side = 1, at=seq(2, length(num.snps), by=2), labels = rev(num.snps[c(TRUE,FALSE)]), tick = FALSE, line=1, cex.axis=0.5)
  if(fit.style=='z') abline(h = -3, lty=2)
  
  spl <- split(m.df, m.df$variable)
  ext.m.df <- do.call(rbind, lapply(spl, function(i){
    switch(fit.style,
           "baf"=i[i$value > quantile(i$value, 0.95),],
           "z"=i[i$value < quantile(i$value, 0.05),])
  }))
  beeswarm(value ~ variable, data = ext.m.df, method = 'swarm', corral = 'gutter',
           cex=0.8, pch = 16, pwcol=clid, add=TRUE)
  dev.off()
}
