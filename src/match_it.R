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
  require(DNAcopy)
  
  ## Set dataset colors
  dataset.cols <- setNames(RColorBrewer::brewer.pal(6, "Dark2"),
                           c("GDSC", "CCLE", "gCSI", "CGP", "Pfizer"))
  
  PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/CCLid/CCLid"
  ref.dat <- CCLid::loadRef(PDIR, 'baf', bin.size=5e5)
  
  
  var.dat <- format.dat$var
  var.dat <- ref.dat$var
  ref.mat <- ref.dat$ref
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
  #x.mat[is.na(x.mat)] <- median(x.mat, na.rm=TRUE)
  
  x.dist <- similarityMatrix(x.mat, 'euclidean')
  D.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=meta.df)
  balanced <- balanceGrps(D.vals)
  
  models <- trainLogit(balanced, predictors=c('baf'))
  x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
  pred <- assemblePredDat(x.vals, known.class=FALSE)
  pred <- mkPredictions(pred, models)
  
  ## Format extra columns on prediction matrix
  pred$Var1 <- as.character(pred$Var1)
  pred$Var2 <- as.character(pred$Var2)
  pred$clA <- as.character(gsub(ds.pattern, "", pred$Var1))
  pred$clB <- as.character(gsub(ds.pattern, "", pred$Var2))
  pred$cvclA <- meta.df[match(pred$clA, meta.df$ID),]$CVCL
  pred$cvclB <- meta.df[match(pred$clB, meta.df$ID),]$CVCL
  pred$g.truth <- with(pred, clA == clB)
  pred$g.truth <- c("NM", "M")[as.integer(pred$g.truth) + 1]
  pred$g.truth <- factor(pred$g.truth, levels=c("M", "NM"))
  pred$baf.p.fit <- factor(pred$baf.p.fit, levels=c("M", "NM"))
  datasets <-c('CCLE', 'GDSC', 'GNE')
  pred$ds <- with(pred, paste(Var1, Var2, sep=";")) %>% 
    gsub("_.*;", ";", .) %>% gsub("_.*", "", .)
  
  ## Split based on datasets
  comb <- cbind(combn(datasets, m = 2),
                matrix(rep(datasets, 2), ncol=3, byrow = TRUE))
  cmb.pred <- apply(comb, 2, function(ds){
    id1 <- paste(ds[1], ds[2], sep=";")
    id2 <- paste(ds[2], ds[1], sep=";")
    idx <- as.logical((pred$ds == id1) + (pred$ds == id2))
    split(pred, f=idx)[['TRUE']]
  })
  names(cmb.pred) <- apply(comb, 2, paste, collapse="-")

  ## Plot the confusion matrix between any two datasets
  pdf(file.path(PDIR, "match_it", "gne-gdsc-ccle_conc.pdf"), height = 12)
  par(mfrow=c(6,3), mar=c(3, 2, 3, 2))
  P.m.nm <- lapply(names(cmb.pred), function(pid){
    p <- cmb.pred[[pid]]
    conf.m <- table(p$baf.p.fit, p$g.truth)
    fourfoldplot(conf.m, space = 0.2)
    # conf.m[2,2] <- sum(p$g.truth=='M')
    conf.m[2,2] <- 600
    fourfoldplot(conf.m, space = 0.2, conf.level = 0, std='ind.max', 
                 main=pid, col=c("#ca0020", "#404040"))
    
    p.m.nm <- CCLid:::splitToMNM(p) #CCLid:::
    p.m.nm$cellosaurus <- CCLid:::checkAgainst(p.m.nm) #CCLid:::
    err.pcl <- CCLid:::genErrBp(p.m.nm) #CCLid:::
    p.m.nm <- cbind(p.m.nm, err.pcl)
    p.m.nm <- p.m.nm[order(err.pcl[,1], err.pcl[,2]),]
    
    print(table(err.pcl))
    barplot(t(table(err.pcl)), horiz=TRUE, las=1, xlim=c(0,80), 
            col=c("FALSE"="#f4a582", "TRUE"="#b2182b"), border=NA)
    return(p.m.nm)
  })
  names(P.m.nm) <- names(cmb.pred)
  dev.off()
  cat("rl ", file.path(PDIR, "match_it", "gne-gdsc-ccle_conc.pdf\n"))

  ## Write otu unknown cases for manual curation
  lapply(names(P.m.nm), function(Pid){
    P <- P.m.nm[[Pid]]
    write.table(P[,-c(6:10)], file=file.path(PDIR, "match_it", paste0("X-", Pid, ".tsv")),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)
  })

  lapply(c("CCLE", "GDSC", "GNE"), function(ds.id){
    ds.cl.ids <- sapply(names(P.m.nm), function(Pid){
      P <- P.m.nm[[Pid]]
      row.idx <- which(P$err=='X' & (!as.logical(P$pcl)))
      cl.ids <- c(P[row.idx,]$Var1, P[row.idx,]$Var2)
      unique(cl.ids[grep(ds.id, cl.ids)])
    })
    cat(paste(unique(sort(unlist(ds.cl.ids))), collapse="\n"))
  })
  
  lapply(names(P.m.nm), function(Pid){
    P <- P.m.nm[[Pid]]
    write.table(P[,-c(6:10)], file=file.path(PDIR, "match_it", paste0("X-", Pid, ".tsv")),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)
  })
}
