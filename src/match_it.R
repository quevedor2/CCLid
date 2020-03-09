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
  # s.range <- sort(sample(1:length(all.vcfs), size=100))
  # f1.scores.by.snps <- lapply(seq(40, 10, by=-2), function(num.snps){
  vcf.f1.scores <- mclapply(all.vcfs[1:400], function(vcf){
    vcf.map <- mapVcf2Affy(file.path(vcf.dir, vcf))
    
    cat(paste0(vcf, "(", match(vcf, all.vcfs), "/", length(all.vcfs), "): "))
    f1.scores <- sapply(num.snps.to.test, function(num.snps){
      cat(paste0(num.snps, "."))
      idx <- sample(1:length(ref.dat$var), size=max(num.snps.to.test)*10, replace = FALSE)
      vcf.map.var <- mapVariantFeat(vcf.map, ref.dat$var[idx])
      vcf.map.var$BAF <- vcf.map.var$BAF[1:num.snps,]
      vcf.map.var$GT <- vcf.map.var$GT[1:num.snps,]
      
      vaf.to.map <- vcf.map.var
      
      ## Overlap the two VCFs to form a single matrix to combine
      ov.idx <- overlapPos(comp = vaf.to.map$BAF,
                           ref=ref.dat$ref, mapping = 'probeset')
      x.mat <- cbind(vaf.to.map$BAF$BAF[ov.idx$comp], 
                     ref.dat$ref[ov.idx$ref,])
      # if(rm.gne){
      #   gne.idx <- c(grep("^GNE_", colnames(x.mat)), grep("^Unk[0-9]", colnames(x.mat)))
      #   x.mat <- x.mat[,-gne.idx]
      # }
      
      ## Calculate distance between samples
      x.dist <- similarityMatrix(x.mat, 'euclidean')
      D.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=meta.df)
      balanced <- balanceGrps(D.vals)
      
      ## Train model
      models <- trainLogit(balanced, predictors=c('baf'))
      x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
      pred <- assemblePredDat(x.vals, known.class=FALSE)
      pred <- mkPredictions(pred, models)
      x.pred <- split(pred, pred$Var2)[[colnames(x.mat)[1]]]
      x.pred <- x.pred[order(x.pred$q),]
      
      rna.idx <- grep(gsub(".snpOut.*", "", vcf), rna.meta.df$EGAF)
      match.idx <- factor(grepl(paste0("_?", rna.meta.df[rna.idx,]$ID, "$"), x.pred$Var1), levels=c(TRUE,FALSE))
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
  save(vcf.f1.scores, file=file.path(PDIR, "match_it", paste0("vcf_f1.", seed, "_", dataset, ".rda")))
  load(file.path(PDIR, "match_it", paste0("vcf_f1.", seed, "_", dataset, ".rda")))
  
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
  scp.path <- "scp quever@192.168.198.99:"
  cat(paste0(scp.path, file.path(vcf.dir, paste0("f1score_", dataset, ".pdf")), ' .\n'))
  
}



##########################
#### Match All RNAseq ####
##########################
## Matches all the RNAseq data from CCLE and GDSC2
expressThis <- function(){
  dataset <- 'CCLE' #GDSC
  ds.pattern = '(CCLE)_|(GDSC)_|(GDSC2)_|(CCLE2)_|(GNE)_|(GNE2)_'
  
  vcf.dir <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs',
                       dataset)
  all.vcfs <- list.files(vcf.dir, pattern="vcf.gz$")
  rna.meta.df <- readinRnaFileMapping()

  rna.identity <- mclapply(all.vcfs, function(vcf, num.snps=200, top.hits=4){
    vcf.map <- mapVcf2Affy(file.path(vcf.dir, vcf))
    
    cat(paste0(vcf, "(", match(vcf, all.vcfs), "/", length(all.vcfs), "): "))
    idx <- sample(1:length(ref.dat$var), size=max(num.snps)*10, replace = FALSE)
    vcf.map.var <- mapVariantFeat(vcf.map, ref.dat$var[idx])
    vcf.map.var$BAF <- vcf.map.var$BAF[1:num.snps,]
    vcf.map.var$GT <- vcf.map.var$GT[1:num.snps,]
    vaf.to.map <- vcf.map.var
    
    ## Overlap the two VCFs to form a single matrix to combine
    ov.idx <- overlapPos(comp = vaf.to.map$BAF,
                         ref=ref.dat$ref, mapping = 'probeset')
    x.mat <- cbind(vaf.to.map$BAF$BAF[ov.idx$comp], 
                   ref.dat$ref[ov.idx$ref,])
    # if(rm.gne){
    #   gne.idx <- c(grep("^GNE_", colnames(x.mat)), grep("^Unk[0-9]", colnames(x.mat)))
    #   x.mat <- x.mat[,-gne.idx]
    # }
    
    ## Calculate distance between samples
    x.dist <- similarityMatrix(x.mat, 'euclidean')
    D.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=meta.df)
    balanced <- balanceGrps(D.vals)
    
    ## Train model
    models <- trainLogit(balanced, predictors=c('baf'))
    x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
    pred <- assemblePredDat(x.vals, known.class=FALSE)
    pred <- mkPredictions(pred, models)
    x.pred <- split(pred, pred$Var2)[[colnames(x.mat)[1]]]
    x.pred <- x.pred[order(x.pred$q),]
    
    if(length(grep("^M$", x.pred$baf.p.fit)) < top.hits){
      return(x.pred[1:top.hits,])
    } else {
      return(x.pred[grep("^M$", x.pred$baf.p.fit), ])
    }
  }, mc.cores = 4)
  names(rna.identity) <- gsub(".snpOut.*", "", all.vcfs)
  err.idx <- sapply(rna.identity, function(i) class(i) == 'try-error')
  if(any(err.idx)) rna.identity <- rna.identity[-which(err.idx)]
  save(rna.identity, file=file.path(PDIR, "match_it", paste0(dataset, "_rnaID.rda")))
  load(file.path(PDIR, "match_it", paste0(dataset, "_rnaID.rda")))
  
  ## Annonate the CVCL IDs
  rna.cvcl <- lapply(names(rna.identity), function(id, dataset){
    print(id)
    rna.id <- rna.identity[[id]]
    file.prefix <- switch(dataset,
                          "GDSC"="EGAF",
                          "CCLE"="SRR")
    if(any(grepl(paste0("^", id, "$"), rna.meta.df[[file.prefix]]))){
      rna.id$Var2 <-  rna.meta.df[grep(paste0("^", id, "$"), rna.meta.df[[file.prefix]]),]$ID
      
      rna.id$clA <- as.character(gsub(ds.pattern, "", rna.id$Var1))
      rna.id$clB <- as.character(gsub(ds.pattern, "", rna.id$Var2))
      rna.id$cvclA <- meta.df[match(rna.id$clA, meta.df$ID),]$CVCL
      rna.id$cvclB <- meta.df[match(rna.id$clB, meta.df$ID),]$CVCL
      rna.id$g.truth <- with(rna.id, clA == clB)
      rna.id[rna.id == 'character(0)'] <- 'CVCL_X482'
      rna.id[is.na(rna.id)] <- 'CVCL_X482'
      if(any(grepl("\\\t$", rna.id$cvclA))) rna.id$cvclA <- gsub("\\\t$", "", rna.id$cvclA)
      if(any(grepl("\\\t$", rna.id$cvclB))) rna.id$cvclB <- gsub("\\\t$", "", rna.id$cvclB)
      rna.id$cellosaurus <- CCLid:::checkAgainst(rna.id) #CCLid:::
      rna.id[order(rna.id$cellosaurus, decreasing=TRUE),]
    }  else {
      rna.id$Var2 <-  id
      rna.id
    }
  }, dataset=dataset)
  names(rna.cvcl)  <- names(rna.identity)
  no.match.idx <- which(sapply(rna.cvcl, function(i) !any(i$g.truth)))
  match.idx <- which(sapply(rna.cvcl, function(i) any(i$g.truth)))
  
  ## Descriptive stats of the number of "matching" and "nonmatching" based on CVCL
  print(paste0("Match: ", length(match.idx),  " / ", length(rna.cvcl)))
  print(paste0("No-match: ", length(no.match.idx),  " / ", length(rna.cvcl)))
  
  ## Concatenate everything and output csv
  rna.m <- plyr::rbind.fill(lapply(rna.cvcl[match.idx], function(i){
    i[which(i$baf.p.fit == 'M'),]
  }))
  rna.nm <- plyr::rbind.fill(rna.cvcl[no.match.idx])
  rna.mnm <- plyr::rbind.fill(rna.m, rna.nm)
  
  table(sapply(rna.cvcl[no.match.idx], function(i) i[1,'baf.p.fit'])) # GDSC-M:14  NM:19 || CCLE-M:13  NM:44 
  table(sapply(rna.cvcl[match.idx], function(i) any(i$g.truth))) # GDSC-M:14  NM:19 || CCLE-M:13  NM:44 
  
  cat(paste0("grep \"", names(rna.cvcl[no.match.idx]), "\" *\n"), "\n")
  
  write.table(rna.mnm, file.path(PDIR, "match_it", paste0(dataset, "_rnaID.csv")),
              sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
  scp.path <- "scp quever@192.168.198.99:"
  cat(paste0(scp.path, file.path(PDIR, "match_it", paste0(dataset, "_rnaID.csv")), ' .\n'))
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
