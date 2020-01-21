bencharkCCLid <- function(bench){
  library(VariantAnnotation)
  library(CCLid)
  
  ## Load in Ref mat file
  PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/CCLid/CCLid"
  analysis <- 'baf'
  ref.mat <- downloadRefCCL("BAF", saveDir = PDIR)
  format.dat <- formatRefMat(name="BAF", ref.mat=ref.mat, 
                             analysis='baf', bin.size=5e5)
  ref.mat <- format.dat$mat
  var.dat <- format.dat$var
  
  ## Load in VCF file
  vcfFile <- '/mnt/work1/users/home2/quever/xfer/A549.sample_id.vcf' ## A549 WES
  vcf.map <- mapVcf2Affy(vcfFile)
}
  
combinedCellLine <- function(){
  new.ids <- assignGrpIDs(ref.mat, meta.df)
  colnames(ref.mat) <- new.ids
  ref.mat.bkup <- ref.mat
  
  a <- 'A549'; b <- 'DU-145'
  cl.A <- grep(paste0("CCLE_", a), colnames(ref.mat))
  cl.B <- grep(paste0("GDSC_", DU-145), colnames(ref.mat))
  sample.mat <- ref.mat[,c(cl.A, cl.B)]
  
  vcf.map.var <- mapVariantFeat(sample.mat, var.dat)
  vcf.to.use <- vcf.map.var[,c(1,2)]
  rownames(vcf.to.use) <- vcf.map.var$Probe_Set_ID

  all.probs <- lapply(seq(0, 1, by=0.05), function(p){
    prop <- c(p, 1-p)
    sample.x <- as.matrix(combineSamples('BAF', vcf.to.use, prop))
    colnames(sample.x) <- "X"

    ov.idx <- overlapPos(comp = sample.x, ref=ref.mat, 
                         mapping = 'probeset')
    x.mat <- cbind(sample.x[ov.idx$comp], 
                   ref.mat[ov.idx$ref,])
    
    
    x.dist <- similarityMatrix(x.mat, 'euclidean')
    as.matrix(c(head(x.dist[order(x.dist[,1,drop=FALSE], decreasing = TRUE),], 10),
                head(x.dist[order(x.dist[,1,drop=FALSE], decreasing = FALSE),], 10)))
    
    
    D.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=meta.df)
    balanced <- balanceGrps(D.vals)
    # plotHist(D.vals$baf)
    
    models <- trainLogit(balanced, predictors=c('baf'))
    
    x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
    pred <- assemblePredDat(x.vals, known.class=FALSE)
    pred <- mkPredictions(pred, models)
    
    x.pred <- split(pred, pred$Var2)[[1]]
    a.prob <- (1 - x.pred[grep(a, x.pred$Var1),]$baf.fit)
    b.prob <- (1- x.pred[grep(gsub("-", ".", b), x.pred$Var1),]$baf.fit)
    names(a.prob) <- paste0(a, "_", c(1:2))
    names(b.prob) <- paste0(b, "_", c(1:2))
    
    c(a.prob, b.prob)
  })
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

snpsCellIdentity <- function(){
  
}