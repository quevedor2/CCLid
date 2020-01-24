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
  new.ids[duplicated(new.ids)] <- gsub("_", "2_",new.ids[duplicated(new.ids)])
  colnames(ref.mat) <- new.ids
  
}

readinRnaFileMapping <- function(){
  meta <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs/fileList1357.txt'
  meta2 <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs/E-MTAB-3983.sdrf.txt'
  
  meta <- read.table(meta, sep="\t", header=F, fill=TRUE)
  meta2 <- read.table(meta2, sep="\t", header=T, fill=TRUE)
  all.meta <- merge(meta, meta2, by.x='V2', by.y='Comment.EGA_SAMPLE.', all=TRUE)
  all.meta <- all.meta[,c('V1','Source.Name', 'V4', 'Comment.EGA_RUN.')]
  
  all.meta$tmp <- gsub("[ -/]", "", all.meta$'Source.Name')
  meta.df$tmp <- gsub("[ -/]", "", meta.df$ID)
  all.meta.df <- merge(all.meta, meta.df, by="tmp", all.x=TRUE)
  colnames(all.meta.df)[1:6] <- c("tmp", "V1", "ID2", "EGAF", "EGAR", "ID")
  return(all.meta.df)
}

###################################
#### Drift between SNP and RNA ####
###################################
## This process is meant to genetic fraction
## of drift found in the RNAseq and the SNP array
## data.  As well as calculate the overlap of
## segments between the two
driftTech <- function(){
  dataset <- 'GDSC'
  vcf.dir <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs',
                       dataset)
  all.vcfs <- list.files(vcf.dir, pattern="vcf.gz$")
  rna.meta.df <- readinRnaFileMapping()
  seed <- 1234
  set.seed(seed)
  
  s.range <- sample(1:length(all.vcfs), size=100)
  # f1.scores.by.snps <- lapply(seq(40, 10, by=-2), function(num.snps){
  vcf.f1.scores <- mclapply(all.vcfs[s.range], function(vcf){
    vcf.map <- mapVcf2Affy(file.path(vcf.dir, vcf))
    vcf.map.var <- mapVariantFeat(vcf.map, var.dat)
    vaf.to.map <- vcf.map.var
    
    ## Overlap the two VCFs to form a single matrix to combine
    ov.idx <- overlapPos(comp = vaf.to.map$BAF,
                         ref=ref.mat, mapping = 'probeset')
    x.mat <- cbind(vaf.to.map$BAF$BAF[ov.idx$comp], 
                   ref.mat[ov.idx$ref,])
    
  })
}
  
############################
#### F1 Score Min.Snps  ####
############################
## This process is meant to calculate the minimum
## number of SNPS needed by calculating the F1 score
## for a series of different number of SNPs
minimumSnps <- function(){
  dataset <- 'GDSC'
  vcf.dir <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs',
                       dataset)
  all.vcfs <- list.files(vcf.dir, pattern="vcf.gz$")
  rna.meta.df <- readinRnaFileMapping()
  num.snps.to.test <- seq(100,10,by=-10)
  seed <- 1234
  
  set.seed(seed)
  s.range <- sort(sample(1:length(all.vcfs), size=100))
  # f1.scores.by.snps <- lapply(seq(40, 10, by=-2), function(num.snps){
  vcf.f1.scores <- mclapply(all.vcfs[s.range], function(vcf){
    vcf.map <- mapVcf2Affy(file.path(vcf.dir, vcf))
    
    cat(paste0(vcf, "(", match(vcf, all.vcfs), "/", length(all.vcfs), "): "))
    f1.scores <- sapply(num.snps.to.test, function(num.snps){
      cat(paste0(num.snps, "."))
      idx <- sort(sample(1:length(var.dat), size=num.snps, replace = FALSE))
      vcf.map.var <- mapVariantFeat(vcf.map, var.dat[idx])
      vaf.to.map <- vcf.map.var
      
      ## Overlap the two VCFs to form a single matrix to combine
      ov.idx <- overlapPos(comp = vaf.to.map$BAF,
                           ref=ref.mat, mapping = 'probeset')
      x.mat <- cbind(vaf.to.map$BAF$BAF[ov.idx$comp], 
                     ref.mat[ov.idx$ref,])
      
      ## Calculate distance between samples
      x.dist <- similarityMatrix(x.mat, 'euclidean')
      D.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=meta.df)
      balanced <- balanceGrps(D.vals)
      
      ## Train model
      models <- trainLogit(balanced, predictors=c('baf'))
      x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
      pred <- assemblePredDat(x.vals, known.class=FALSE)
      pred <- mkPredictions(pred, models)
      x.pred <- split(pred, pred$Var2)[[1]]
      
      rna.idx <- grep(gsub(".snpOut.*", "", vcf), rna.meta.df$EGAF)
      match.idx <- grepl(paste0("_", rna.meta.df[rna.idx,]$ID, "$"), x.pred$Var1)
      c.tbl <- sapply(split(x.pred, match.idx), function(i) table(i$baf.p.fit))
      
      if(all(dim(c.tbl) == c(2,2))){
        c.tbl <- c.tbl[c("M", "NM"), c('TRUE', 'FALSE')]
        
        precision <- c.tbl[1,1] / (c.tbl[1,1] + c.tbl[1,2])
        recall <- c.tbl[1,1] / (c.tbl[1,1] + c.tbl[2,1])
        f1.score <- 2 * ((precision * recall) / (precision + recall))
        if(is.nan(f1.score)) f1.score <- 0
      } else {
        f1.score <- NULL
      }
      return(f1.score)
    })
    cat("\n")
    
    gc()
    return(f1.scores)
  }, mc.cores = 4)
  save(vcf.f1.scores, file=file.path(vcf.dir, paste0("vcf_f1.", seed, "_", dataset, ".rda")))
  
  pdf(file.path(vcf.dir, paste0("f1score_", dataset, ".pdf")), height = 4, width=5)
  f1.mat <- do.call(rbind, lapply(vcf.f1.scores, unlist))
  colnames(f1.mat) <- num.snps.to.test
  dataset <- 'GDSC'
  boxplot(f1.mat, col="#0000ff22", las=1, outline=FALSE, 
          ylab='F1-Score', xlab='Num. of SNPs', main=dataset)
  rm.idx <- apply(f1.mat, 2, function(i) i >= median(i))
  f1.mat[rm.idx] <- NA
  colnames(f1.mat) <- rev(colnames(f1.mat))
  melt.f1 <- reshape::melt(f1.mat)
  
  beeswarm::beeswarm(value ~ X2, data = melt.f1, method = 'swarm', corral = 'gutter',
           cex=0.8, pch = 16, add=TRUE, col=scales::alpha("black", 0.6))
  dev.off()
}
