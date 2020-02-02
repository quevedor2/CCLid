readinRnaFileMapping <- function(){
  meta <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/data/GDSC/fileList1357.txt'
  meta.gdsc <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/data/GDSC/E-MTAB-3983.sdrf.txt'
  meta.ccle <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/data/CCLE/CCLE_meta.txt'
  pattern="[-\\(\\)\\.\\,\\_\\/ ]"
  
  meta <- read.table(meta, sep="\t", header=F, fill=TRUE)
  meta.gdsc <- read.table(meta.gdsc, sep="\t", header=T, fill=TRUE)
  meta.ccle <- read.table(meta.ccle, sep=",", header=T, fill=TRUE)
  
  meta.gdsc$simpleid = toupper(gsub(pattern, "", meta.gdsc$Source.Name))
  meta.ccle$simpleid = toupper(gsub(pattern, "", meta.ccle$Cell_Line))
  ov = sort(intersect(meta.gdsc$simpleid, meta.ccle$simpleid))  ## 61 from non simple.id, 79 simple
  gdsc.o = sort(setdiff(meta.gdsc$simpleid, meta.ccle$simpleid))
  ccle.o = sort(setdiff(meta.ccle$simpleid, meta.gdsc$simpleid))
  
  # Merge by EGAF(meta) to EGAR (meta.gdsc) and cell-name by EGAN id
  all.meta <- merge(meta, meta.gdsc, by.x='V2', by.y='Comment.EGA_SAMPLE.', all=TRUE)
  all.meta <- merge(all.meta, meta.ccle, by='simpleid', all=TRUE)
  all.meta <- all.meta[,c('V1','Source.Name', 'V4', 'Comment.EGA_RUN.', 'Run', 
                          'Cell_Line', 'simpleid')]
  
  meta.df$simpleid <- gsub(pattern, "", meta.df$ID)
  meta.df[grep("^T-T$", meta.df$ID),]$simpleid <- 'T-T'
  all.meta.df <- merge(all.meta, meta.df, by="simpleid", all.x=TRUE)
  
  colnames(all.meta.df)[1:8] <- c("tmp", "V1", "GDSC_ID", "EGAF", "EGAR", "SRR", "CCLE_ID", "ID")
  return(all.meta.df)
}


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
  format.dat <- formatRefMat(name="BAF", ref.mat=ref.mat, saveDir=PDIR, 
                             analysis='baf', bin.size=5e5)
  ref.mat <- format.dat$mat
  var.dat <- format.dat$var
  rm(format.dat)
  
  new.ids <- assignGrpIDs(ref.mat, meta.df)
  new.ids[duplicated(new.ids)] <- gsub("_", "2_",new.ids[duplicated(new.ids)])
  colnames(ref.mat) <- new.ids
}

######################################################
#### Genotype concordance between every cell line ####
######################################################
## This process is meant to compare the genotypes using
## BAF between every cell line in all datasets used
snpsCellIdentity <- function(){
  vcf.map.var <- mapVariantFeat(ref.mat, var.dat)
  x.mat <- as.matrix(vcf.map.var)
  storage.mode(x.mat) <- 'numeric'
  ds.pattern = '(CCLE)_|(GDSC)_|(GDSC2)_|(CCLE2)_|(GNE)_|(GNE2)_'
  rm.idx <- grep(ds.pattern, colnames(x.mat), invert = TRUE)
  x.mat <- x.mat[,-rm.idx]
  
  boxplot(apply(x.mat, 1, function(i) { sum(is.na(i))}))
  x.dist <- similarityMatrix(x.mat, 'euclidean')
  D.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=meta.df)
  balanced <- balanceGrps(D.vals)
  
  models <- trainLogit(balanced, predictors=c('baf'))
  x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
  pred <- assemblePredDat(x.vals, known.class=FALSE)
  pred <- mkPredictions(pred, models)
  
  
  pred$g.truth <- gsub(ds.pattern, "", pred$Var1) == gsub(ds.pattern, "", pred$Var2)
  pred$g.truth <- c("NM", "M")[as.integer(pred$g.truth) + 1]
  conf.m <- table(pred$baf.p.fit, pred$g.truth)
  fourfoldplot(conf.m, space = 0.5)
  
  ## Check for clear isogenic variants based on annotation
  p.FN$strdist <- apply(p.FN, 1, function(i) {
    a=toupper(gsub("[^a-zA-Z0-9]", "", gsub(ds.pattern, "", i['Var1'])))
    b=toupper(gsub("[^a-zA-Z0-9]", "", gsub(ds.pattern, "", i['Var2'])))
    if(nchar(a) > nchar(b)) grepl(b,a) else grepl(a,b)
  })
  
  pred2 <- pred
  pred2$Var1 <- gsub("GDSC", "GNE", pred2$Var1)
  pred2$Var2 <- gsub("GDSC", "GNE", pred2$Var2)
  p.FN <- .getFN(pred, ds.pattern=ds.pattern)
  p2.FN <- .getFN(pred2, ds.pattern=ds.pattern)
  
  pG.FN <- list("G1"=p.FN, "G2"=p2.FN)
  unique(as.character(sapply(pG.FN, function(i) i$pair)))
  
  .simplePair <- function(p, ds.pattern){
    CLid <- data.frame("A"=gsub(ds.pattern, "", p$Var1),
                       "B"=gsub(ds.pattern, "", p$Var2))
    apply(CLid, 1, function(i) paste(sort(i), collapse="_"))
  }
  
  .getFN <- function(p, ...){
    FN <- intersect(which(p$baf.p.fit=='M'), which(p$g.truth=='NM'))
    FN.df <- p[FN,]
    FN.df <- FN.df[order(FN.df$baf.p.fit),]
    FN.df$pair <- .simplePair(FN.df, ...)
    return(FN.df)
  }
  
  FN <- intersect(which(pred$baf.p.fit=='M'), which(pred$g.truth=='NM'))
  FN.df <- pred[FN,]
  FN.df <- FN.df[order(FN.df$baf.p.fit),]
  unmapped <- unique(c(grep(ds.pattern, FN.df$Var1, invert = TRUE), grep(ds.pattern, FN.df$Var2, invert = TRUE)))
  FN.df <- FN.df[-unmapped,]
  head(FN.df); tail(FN.df)
  head(pred[FN,][order(pred[FN,]$baf.fit),1:5], 100)
  
  
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
