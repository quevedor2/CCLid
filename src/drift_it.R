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

isolateL2Rpsets <- function(){
  setwd("/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn/log2r_bins")
  for(i in c("CCLE", "GDSC")){
    require(Biobase)
    pset=readRDS(file.path("..", paste0(i, "_CN.bins.RDS")))
    env <- new.env()
    env[['exprs']] <- assayData(pset)$Log2Ratio
    assayData(pset) <- env
    saveRDS(pset, file=paste0(i, "_CN.l2r_bins.RDS"))
  }
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
}


####################################################
#### Concordance between BAF-drift and CN-drift ####
####################################################
## This process is meant to compare the drift called between
## any two matching cell line IDs using both BAF-drift
## from the CCLid package, and the difference
## in ASCAT ASCN data.
driftConcordance <- function(){
  dataset <- 'GDSC'
  alt.ds <- 'CCLE'
  
  ## Find variant features and isolate for cell lines shared in datasets
  ref.mat.var <- mapVariantFeat(ref.dat$ref, ref.dat$var)
  cl.idx <- CCLid::findCclPairs(meta.df, ref.mat.var, ds=c('CCLE', 'GDSC'))
  m.cls.idx <- cl.idx[sapply(cl.idx, function(i) length(i) >= 2)]
  
  ## Get drift distance between CL pairs using BAF
  baf.drifts <- lapply(m.cls.idx, CCLid::getBafDrifts, x.mat=ref.mat.var, 
                       ref.ds=dataset, alt.ds=alt.ds)
  save(baf.drifts, file=file.path(PDIR, "drift_it", 
                             paste0(dataset, "-", alt.ds, "_baf_drift.rda")))
  
  ## Load in CN bins
  cn.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn/50kb_bins'
  cn.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn'
  cn.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn/log2r_bins'
  bins.file <- list.files(cn.dir, pattern="bins", include.dirs = FALSE)
  bins <- lapply(bins.file, function(b) readRDS(file.path(cn.dir, b)))
  names(bins) <- gsub("_.*", "",  bins.file)
  
  cn.drift=CCLid::getCNDrifts(ref.l2r=assayData(cn.bins[[dataset]])$exprs,
                              alt.l2r=assayData(cn.bins[[alt.ds]])$exprs,
                              cell.ids=names(baf.drifts), segmenter='PCF')
  save(cn.drifts, file=file.path(PDIR, "drift_it", 
                                paste0(dataset, "-", alt.ds, "_cn_drift.rda")))
  
  ## Find overlap between GRanges of BAF to CN drift
  
  ## Extract and combine BAF drifts and CN drifts
  cn.drift.frac <- reshape::melt(sapply(cn.drifts$frac, function(i) i[[3]]))
  cn.drift.frac$L1 <- rownames(cn.drift.frac)
  baf.drift.frac <- sapply(baf.drifts, function(i) i$frac)
  null.idx <- which(sapply(baf.drift.frac, is.null))
  baf.drift.frac[null.idx] <- NA
  baf.drift.frac <- reshape::melt(baf.drift.frac)
  
  ## Merge the BAF and CN data
  cn.baf.drift <- merge(cn.drift.frac, baf.drift.frac, by="L1", all=TRUE)
  rownames(cn.baf.drift) <- cn.baf.drift$L1
  cn.baf.drift <- cn.baf.drift[,-1]
  colnames(cn.baf.drift) <- c("CN", "BAF")
  
  ## Plot the correlation
  plot(cn.baf.drift, xlim=c(0,1), ylim=c(0,1), main="Fraction of genome drifted")
  r <- round(cor(cn.baf.drift, use = "complete.obs")[1,2],2)
  abline(coef=c(0,1), col="grey", lty=2)
  
  ## Extract and combine BAF drifts and CN drifts
  
  
  
  
  
  [[1]]$
  names(cn.drift)
  cn.drift.frac <- unlist(sapply(cn.drift, CCLid:::.getDrift, idx=1))
  baf.drift.frac <- unlist(sapply(baf.drifts[names(cn.drift)], function(i) i$fraca))
  df <- data.frame("cn"=cn.drift.frac,
                   "baf"=baf.drift.frac)
  
  CCLid:::plot.CCLid(cn.drift[['42-MG-BA']]$cna.obj)
  CCLid:::plot.CCLid(bdf$cna.obj[[1]])
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
  alt.ds <- 'CCLE'
  
  vcf.dir <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs',
                       dataset)
  all.vcfs <- list.files(vcf.dir, pattern="vcf.gz$")
  rna.meta.df <- readinRnaFileMapping()
  names(all.vcfs) <- gsub(".snpOut.*", "", all.vcfs)
  
  ## Compare every VCF to the entire ref matrix to calculate BAF drift
  vcf.drift <- mclapply(all.vcfs, function(vcf){
    cat(paste0(vcf, "(", match(vcf, all.vcfs), "/", length(all.vcfs), "): "))
    vcf.map <- mapVcf2Affy(file.path(vcf.dir, vcf))
    vcf.map.var <- mapVariantFeat(vcf.map, var.dat)
    vaf.to.map <- vcf.map.var
    
    ## Overlap the two VCFs to form a single matrix to combine
    ov.idx <- overlapPos(comp = vaf.to.map$BAF,
                         ref=ref.mat, mapping = 'probeset')
    x.mat <- cbind(vaf.to.map$BAF$BAF[ov.idx$comp], 
                   ref.mat[ov.idx$ref,])
    
    ## Identify matching cell line data to RNAseq
    rna.idx <- grep(gsub(".snpOut.*", "", vcf), rna.meta.df$EGAF)
    match.idx <- grep(paste0("_", rna.meta.df[rna.idx,]$ID, "$"), colnames(x.mat))
    colnames(x.mat)[1] <- paste0("RNA_", rna.meta.df[rna.idx, 'ID'])
    
    ## Calculate drift of Cell line with RNAseq with external control
    print(match.idx)
    if(length(match.idx) > 0){
      x.drift <- bafDrift(sample.mat=x.mat[,c(1, match.idx)])

      ## Isolate siginificant different regions
      sig.diff.gr <- lapply(x.drift$cna.obj, sigDiffBaf)
      frac.cnt <- x.drift$frac
      
      summ <- list("frac"=frac.cnt,
                   "sig"=sig.diff.gr)
    } else {
      summ=NULL
    }
    
    
    return(summ)
  }, mc.cores = 8)
  save(vcf.drift, file=paste0("~/", dataset, "_vcf_drift.rda"))

  ## Identify drift pairs that are NULL and remove from future analysis
  seg.sig <- lapply(vcf.drift, function(i) if(length(i$sig) == 2) i$sig else NULL)
  null.idx <- which(sapply(seg.sig, is.null))

  ## Intersect significant (z > 3) drifted regions and calculate estimates
  ## of how much intersect there is between SNP and RNA data
  all.drift.ovs <- lapply(seg.sig[-null.idx], function(seg){
    drift.ovs <- driftOverlap(seg, ref.ds=dataset, alt.ds=alt.ds)
  })
  null.idx2 <- sapply(all.drift.ovs, is.null)
  if(any(null.idx2)) all.drift.ovs <- all.drift.ovs[-which(null.idx2)]

  # Reduce drift fraction data to a matrix and order
  idx.drift <- lapply(seq_along(all.drift.ovs[[1]]), function(idx){
    drift.mat <- do.call(cbind, lapply(all.drift.ovs, function(i) i[[idx]]))
    colnames(drift.mat) <- names(all.drift.ovs)
    drift.mat
  })
  tmp.m <- idx.drift[[1]]
  ord <- order(tmp.m[3,], tmp.m[1,], tmp.m[2,])
  
  ## Aggregate the "z > 3" drift "fraction of the genome" data into a singular matrix
  rna.drift <- do.call(gtools::smartbind, lapply(vcf.drift, CCLid:::.getDrift, idx=1)) #RNA - GDSC/CCLE
  cl.drift <- do.call(gtools::smartbind, lapply(vcf.drift, CCLid:::.getDrift, idx=2)) # GDSC - CCLE
  colnames(rna.drift) <- paste0("RNA_", colnames(rna.drift))
  colnames(cl.drift) <- paste0(dataset, "_", colnames(cl.drift))
  rna.drift$ID <- rownames(rna.drift)
  cl.drift$ID <- rownames(cl.drift)
  rna.cl.drift <- merge(rna.drift, cl.drift, by="ID", all.y=TRUE)
  rcl.sub.drift <- rna.cl.drift[match(names(all.drift.ovs), rna.cl.drift$ID),]
  
  ## Assign the colours and labels
  all.row.ids <- unique(as.character(sapply(idx.drift, rownames)))
  cols <- setNames(RColorBrewer::brewer.pal(length(all.row.ids), "Set2"),
                   all.row.ids)
  cols['intersect'] <- '#ffffcc'
  cols['no_drift'] <- 'grey'

  #### Visualization ####
  ## Plot amount of intra- and inter-institutional drift
  pdf(file.path(vcf.dir, paste0("drift_", dataset, "-", alt.ds, ".pdf")), 
      height = 3, width=5)
  par(mfrow=c(3,1), mar=c(0.5, 4.1, 0.5, 2.1))
  ids <- combn(c("RNA", dataset, alt.ds), 2)
  apply(ids, 2, function(id){
    barplot(rcl.sub.drift[ord, paste(id, collapse="_")], ylab=paste(id, collapse="/"),
            las=1, ylim=c(0,1), col=cols[paste(id, collapse="/")], xaxt='n', xlab='')
  })
  
  ## Plot individual cell lines
  cl.idx <- c(which.max(idx.drift[[1]][ord,]['intersect',]), 
              41, which.max(rcl.sub.drift[ord,]$RNA_GDSC))
  par(mfrow=c(3,1), mar=c(0.5, 4.1, 0.5, 2.1))
  sapply(colnames(idx.drift[[1]][,ord])[cl.idx], function(cl.ids){
    #rna.meta.df[sapply(colnames(idx.drift[[1]][,ord])[cl.idx], grep, rna.meta.df$EGAF),]$ID
    multiDriftPlot(seg.sig[-null.idx][[cl.ids]], ref.ds=dataset, alt.ds=alt.ds)
  })
  
  ## Plot the drift overlap matrix
  par(mfrow=c(3,1), mar=c(2, 4.1, 2, 2.1))
  lapply(idx.drift, function(m){
    bp <- barplot(as.matrix(m[,ord]), xaxt='n', col=cols[rownames(m)], las=1)
    
    par(xpd=TRUE)
    if(all(rownames(m) == rownames(idx.drift[[1]]))){
      idx <- match(colnames(m[,ord]), rna.meta.df$EGAF)
      lbl <- rep(NA, ncol(m))
      lbl[cl.idx] <- rna.meta.df[idx, ][cl.idx,]$ID
      axis(side = 1, at=bp, labels=lbl, las=2, cex.axis=0.8)
    } else {
      legend(x = 1, y = 0, fill=cols[rownames(m)], 
             legend=rownames(m), horiz=TRUE, box.lwd = 0)
    }
    par(xpd=FALSE)
  })

  ## Plot concordance in drift estimates
  col.idx <- grep(paste0(alt.ds, "$"), colnames(rna.cl.drift))
  comp.drift <- rna.cl.drift[,col.idx,drop=FALSE]
  par(mfrow=c(2,1), mar=c(2, 4.1, 2, 2.1))
  plot(comp.drift, pch=16, col=scales::alpha('black', 0.7), 
       xlim=c(0,1), ylim=c(0,1),
       xlab=paste0(dataset, " (", gsub("_.*", "", colnames(comp.drift)[1]), ") - ", alt.ds, " (SNP)"),
       ylab=paste0(dataset, " (", gsub("_.*", "", colnames(comp.drift)[2]), ") - ", alt.ds, " (SNP)"),
       main='Genetic drift between cell lines (Fraction of genome)',
       cex=0.6, las=1)
  r <- round(cor(comp.drift, use="pairwise", method="pearson")[1,2], 2)
  text(c(1,0.9), labels = paste0("r = ", r, " (pearson)"), adj=1, cex=0.6)
  abline(coef = c(0,1), col="grey", lty=2)

  dev.off()
  print(file.path(vcf.dir, paste0("drift_", dataset, "-", alt.ds, ".pdf")))

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
