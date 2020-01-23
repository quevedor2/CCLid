bencharkCCLid <- function(bench){
  library(VariantAnnotation)
  library(CCLid)
  
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
  
combinedCellLine <- function(){
  ref.mat.bkup <- ref.mat
  new.ids <- assignGrpIDs(ref.mat.bkup, meta.df)
  colnames(ref.mat) <- new.ids
  
  a <- 'A549'; b <- 'DU-145'
  cl.A <- grep(paste0("CCLE_", a), colnames(ref.mat))
  cl.B <- grep(paste0("GDSC_", b), colnames(ref.mat))
  sample.mat <- ref.mat[,c(cl.A, cl.B)]
  
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
    
    ## Logit model based on euclidean distance between probesets
    analyze.probesets=FALSE
    if(analyze.probesets){
      snpD.vals <- lapply(list("baf"=x.mat), splitSnpDist, meta.df=meta.df)
      balanced <- snpD.vals$baf
      balanced$id <- factor(balanced$id)
      probeset.model <- trainLogit(balanced, predictors=c('.'))
      
      pred.snp <- t(sapply(c(1:ncol(x.mat)), function(c.idx){
        abs(rowDiffs(as.matrix(x.mat[,c(1,c.idx)])))
      }))
      colnames(pred.snp) <- rownames(x.mat)
      rownames(pred.snp) <- colnames(x.mat)
      x.pred.snp <- predict(probeset.model[[1]], newdata = as.data.frame(pred.snp), type = "response")
      
      head(as.matrix(sort(x.pred.snp)), 30)
      as.matrix(sort(x.pred.snp[c('CCLE_DU-145', 'GDSC_DU-145',
                                  'CCLE_A549', 'GDSC_A549')]))
      head(as.matrix(sort(rowSums(pred.snp))), 10)
      lapply(split(rowSums(snpD.vals$baf[,-ncol(snpD.vals$baf)]), snpD.vals$baf$id), summary)
      
      pred <- mkPredictions(pred.snp, probeset.model)
    }
    ## Logit model based on euclidean distance between samples
    x.dist <- similarityMatrix(x.mat, 'euclidean')
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
    truth.cl <- lapply(setNames(c(a,b), c(a,b)), function(i){
      i.idx <- grep(i, x.pred$Var1)
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
    # pred.cl[['OVSAHO']] <- x[3,,drop=FALSE]
    # rownames(pred.cl[['OVSAHO']]) <- "CCLE_OVSAHO"
    A <- sample.x[ov.idx$comp,,drop=FALSE]
    A2 <- A + runif(n=nrow(A), min = -0.2, max = 0.2)
    M <- sapply(pred.cl, function(p.cl){
        combineSamples("BAF", ref.mat[ov.idx$ref, rownames(p.cl), drop=FALSE], 
                     rep(1/nrow(p.cl), nrow(p.cl)))
    })
    M.mse <- checkMse(A2, M)
    
    # Decomposition with complete data
    #M1.mse <- .checkMse(A, M)
    M0 <- NNLM::nnmf(A2, k = 0, check.k = FALSE, init=list(W0 = M));
    M1 <- NNLM::nnmf(A2, k = 1, check.k = FALSE, init=list(W0 = M));
    M1.deconv <- as.matrix(M1$H[,1] / colSums(M1$H)) # [0.9, 0.09, 0.01]
    M0.deconv <- as.matrix(M0$H[,1] / colSums(M0$H)) # [0.9, 0.09, 0.01]
    
    return(list("truth"=truth.cl,
                "pred"=pred.cl,
                "nmf"=M0.deconv))
  })
  lapply(all.deconv, function(i) i$pred)
  lapply(all.deconv, function(i) i$nmf)
  
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

combinedCellLine <- function(){
  ref.mat.bkup <- ref.mat
  new.ids <- assignGrpIDs(ref.mat, meta.df)
  colnames(ref.mat) <- new.ids
  
  
  vcf.dir <- '/mnt/work1/users/home2/quever/xfer/kelsie'
  data.type <- 'wes'
  
  vcf.files <- switch(data.type,
                      "wes"=c('DU-145.sample_id.vcf', 'DU-R600.vcf'), #A549.sample_id.vcf [not A549]
                      "rna"=c('PAR_ATTACTCG.vcf', 'R400_AGCGATAG.vcf'))
   
  ## Load in two VCFs (A549 and DU-145) from a given technology
  vcfs <- lapply(vcf.files, function(v){
    vcf.map <- mapVcf2Affy(file.path(vcf.dir, data.type, v))
    vcf.map.var <- mapVariantFeat(vcf.map, var.dat)
    vcf.map.var
  })
  
  ## Overlap the two VCFs to form a single matrix to combine
  ov.idx <- overlapPos(comp = vcfs[[1]]$BAF,
                       ref=vcfs[[1]]$BAF, mapping = 'probeset')
  vcfs.c <- cbind(vcfs[[1]]$BAF$BAF[ov.idx$comp], 
                  vcfs[[2]]$BAF$BAF[ov.idx$ref])
  colnames(vcfs.c) <- gsub(".vcf", "", vcf.files)
  rownames(vcfs.c) <- vcfs[[1]]$BAF$Probe_Set_ID[ov.idx$comp]
  
  
  
  ## Purpose: Combine the two samples at different proportions
  ## See whether we can detect the two samples and deconvolute them
  lapply(seq(0, 1, by=0.05), function(p){
    prop <- c(p, 1-p)
    sample.x <- as.matrix(combineSamples('BAF', vcfs.c, prop))
    colnames(sample.x) <- "X"
    
    ov.idx <- overlapPos(comp = sample.x, ref=ref.mat, 
                         mapping = 'probeset')
    x.mat <- cbind(sample.x[ov.idx$comp], 
                   ref.mat[ov.idx$ref,])
    
    
    x.dist <- similarityMatrix(x.mat, 'cor')[,1,drop=FALSE]
    as.matrix(c(head(x.dist[order(x.dist, decreasing = TRUE),], 10),
                head(x.dist[order(x.dist, decreasing = FALSE),], 10)))
    
    
    # x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
    # pred <- assemblePredDat(x.vals, known.class=FALSE)
    # pred <- mkPredictions(pred, models)
    # head(pred[order(pred$baf.fit),], n=15)
    # ggplot(pred,aes(x=baf,y=baf.fit)) + stat_binhex()
  })
  
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
  frac.p <- as.data.frame(do.call(rbind, frac.score))
  frac.p$frac.snps <- r
  frac.p$numb.of.snps <- r * length(var.dat)
  
  if(!exists("frac.l")) frac.l <- list()
  frac.l[[data.type]] <- frac.p
  
  pdf("~/minSnps.pdf")
  lapply(names(frac.l), function(data.type){
    plot(0, type="n", xlim=c(min(r), max(r)), ylim=c(0,1), xaxt='n', xlab="Frac. of SNPs", 
         ylab="P", las=1, main=paste0(data.type, ": classify cell identity"))
    
    cols <- setNames(c("#91bfdb", "#4575b4"),
                     colnames(frac.l[[data.type]])[1:2])
    axis(side = 1, at=frac.l[[data.type]]$frac.snps, 
         labels=frac.l[[data.type]]$frac.snps,line=0, tick = TRUE)
    axis(side = 1, at=frac.l[[data.type]]$frac.snps, 
         labels=floor(frac.l[[data.type]]$numb.of.snps), line=1, tick = FALSE)
    for(i in colnames(frac.l[[data.type]])[1:2]){
      points(x = frac.l[[data.type]]$frac.snps, y=frac.l[[data.type]][,i], pch=16, col=cols[i])
      lines(x = frac.l[[data.type]]$frac.snps, y=frac.l[[data.type]][,i], col=cols[i])
    }
    
    legend(x=0.001, y=0.9, fill=cols, legend=names(cols), box.lwd = 0)
  })
  dev.off()
}

snpsCellIdentity <- function(){
  ref.mat.bkup <- ref.mat
  new.ids <- assignGrpIDs(ref.mat, meta.df)
  colnames(ref.mat) <- new.ids
  
  vcf.dir <- '/mnt/work1/users/home2/quever/xfer/kelsie'
  data.type <- 'rna'
  
  vcf.file <- switch(data.type,
                      "wes"=setNames(c('DU-145.sample_id.vcf'), c("DU-145")), #A549.sample_id.vcf [not A549]
                      "rna"=setNames(c('PAR_ATTACTCG.vcf'), c("A549"))) #A549
  
  ## Load in two VCFs (A549 and DU-145) from a given technology
  vcf.map <- mapVcf2Affy(file.path(vcf.dir, data.type, vcf.file))
  
  r <- seq(0.05, 0.001, by=-0.001)
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
    
    x.pred <- split(pred, pred$Var2)[[1]]
    #tail(x.pred[order(x.pred$baf, decreasing = TRUE),])
    cl.id <- names(vcf.file)
    cl.prob <- (1- x.pred[grep(gsub("-", ".", cl.id), x.pred$Var1),]$baf.fit)
    cl.prob <- setNames(cl.prob, as.character(x.pred[grep(gsub("-", ".", cl.id), x.pred$Var1),]$Var1))
    
    cl.prob
  })
  frac.p <- as.data.frame(do.call(rbind, frac.score))
  frac.p$frac.snps <- r
  frac.p$numb.of.snps <- r * length(var.dat)
  
  if(!exists("frac.l")) frac.l <- list()
  frac.l[[data.type]] <- frac.p
  
  pdf("~/minSnps.pdf")
  lapply(names(frac.l), function(data.type){
    plot(0, type="n", xlim=c(min(r), max(r)), ylim=c(0,1), xaxt='n', xlab="Frac. of SNPs", 
         ylab="P", las=1, main=paste0(data.type, ": classify cell identity"))
    
    cols <- setNames(c("#91bfdb", "#4575b4"),
                     colnames(frac.l[[data.type]])[1:2])
    axis(side = 1, at=frac.l[[data.type]]$frac.snps, 
         labels=frac.l[[data.type]]$frac.snps,line=0, tick = TRUE)
    axis(side = 1, at=frac.l[[data.type]]$frac.snps, 
         labels=floor(frac.l[[data.type]]$numb.of.snps), line=1, tick = FALSE)
    for(i in colnames(frac.l[[data.type]])[1:2]){
      points(x = frac.l[[data.type]]$frac.snps, y=frac.l[[data.type]][,i], pch=16, col=cols[i])
      lines(x = frac.l[[data.type]]$frac.snps, y=frac.l[[data.type]][,i], col=cols[i])
    }
    
    legend(x=0.001, y=0.9, fill=cols, legend=names(cols), box.lwd = 0)
  })
  dev.off()
}