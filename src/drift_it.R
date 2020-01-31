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
  format.dat <- formatRefMat(name="BAF", ref.mat=ref.mat, 
                             analysis='baf', bin.size=5e5)
  ref.mat <- format.dat$mat
  var.dat <- format.dat$var
  rm(format.dat)
  
  new.ids <- assignGrpIDs(ref.mat, meta.df)
  new.ids[duplicated(new.ids)] <- gsub("_", "2_",new.ids[duplicated(new.ids)])
  colnames(ref.mat) <- new.ids
  
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
  ref.mat.var <- mapVariantFeat(ref.mat, var.dat)
  cl.idx <- findCclPairs(meta.df, ref.mat.var)
  m.cls.idx <- cl.idx[sapply(cl.idx, function(i) length(i) >= 2)]
  
  ## Get drift distance between CL pairs
  baf.drifts <- lapply(m.cls.idx, function(cl.pairs){
    ref.idx <- grep(paste0(dataset, "_"), colnames(ref.mat.var)[cl.pairs])
    alt.idx <- grep(paste0(alt.ds, "_"), colnames(ref.mat.var)[cl.pairs])
    all.idx <- c(ref.idx, alt.idx)
    
    if(length(all.idx) == 2){
      bdf <- bafDrift(ref.mat.var[,cl.pairs[all.idx]])
      #CCLid:::plot.CCLid(bdf$cna.obj[[1]])
      drift.score <- list("sig.gr"=sigDiffBaf(bdf$cna.obj[[1]]),
                          "frac"=bdf$frac[[1]][3,])
    } else {
      drift.score <- list("sig.gr"=NULL, "frac"=NULL)
    }

    return(drift.score)
  })
  save(baf.drifts, "~/baf_drift.rda")
  # cl.pairs <- m.cls.idx[[1]]
  
  ## Load in CN bins
  cn.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn/50kb_bins'
  cn.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn'
  bins.file <- list.files(cn.dir, pattern="bins", include.dirs = FALSE)[-1]
  bins <- lapply(bins.file, function(b) readRDS(file.path(cn.dir, b)))
  names(bins) <- gsub("_.*", "",  bins.file)
  
  # "exprs"  "nAraw"  "nBraw"  "nMajor" "nMinor" "TCN"   
  alt.bin.ids <- assignGrpIDs(assayData(bins[[alt.ds]])$exprs, meta.df)
  ref.bin.ids <- assignGrpIDs(assayData(bins[[dataset]])$exprs, meta.df)
  alt.ref.idx <- data.frame("id"=as.character(names(baf.drifts)),
                            "ref"=as.integer(sapply(paste0("_", names(baf.drifts), "$"), 
                                                    grep, x=ref.bin.ids)),
                            "alt"=as.integer(sapply(paste0("_", names(baf.drifts), "$"), 
                                                    grep, x=alt.bin.ids)))
  na.idx <- apply(alt.ref.idx, 1, function(i)  any(is.na(i)))
  if(any(na.idx)) alt.ref.idx <- alt.ref.idx[-which(na.idx),]
  cn.drift <- apply(alt.ref.idx, 1, function(ar.i){
    ar.i = unlist(alt.ref.idx[5,])
    ref.id = ref.bin.ids[ar.i['ref']]
    alt.id = alt.bin.ids[ar.i['alt']]
    for( id in list(ref.id, alt.id)){
      did = gsub("_.*", "", id)
      cat(paste("xmor", 
                  file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines',
                            did, 'eacon', names(id), 'ASCAT/L2R', paste0(names(id), '.SEG.ASCAT.png\n'))))
    }
    require(DNAcopy)
    ar.i <- setNames(as.integer(ar.i), names(ar.i))
    print(paste(ar.i, collapse="_"))
    ra.i <- cbind(assayData(bins[[dataset]])$exprs[,ar.i['ref'], drop=FALSE],
                  assayData(bins[[alt.ds]])$exprs[,ar.i['alt'], drop=FALSE])
    idx <- sample(1:nrow(ra.i), size=1000, replace=FALSE)
    # ra.i = scale(ra.i, scale=FALSE)
    # #scaling <- apply(ra.i, 2, norm_vec)
    # scaling <- sapply(1:100, function(i) apply(ra.i[sample(nrow(ra.i), 30000),], 2, norm_vec))
    # scaling <- apply(scaling, 2, function(i) i[1]/i[2])
    # scaling <- scaling[which.min(abs(scaling-1))]
    # ra.i[,1] <- ra.i[,1] / scaling # (scaling[1]/scaling[2])
    ra.i <-as.data.frame(ra.i)
    colnames(ra.i) <-c('ref', 'alt')
    
    plot(ra.i[idx,], xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5))
    abline(coef=c(0,1))
    
    #ra.lm <- lm(alt ~ ref, data=ra.i)
    
    
    #ra.i <- preprocessCore::normalize.quantiles(data.matrix(ra.i))
    D = (ra.i[,1] - ra.i[,2])
    # D[is.na(D)] <- median(D,na.rm=TRUE)
    # D <- scale(D, scale = FALSE)
    # D <- ra.lm$residuals

    CNAo <- with(featureData(bins[[alt.ds]])@data, #[names(ra.lm$residuals),],
                 CNA(genomdat=D, 
                     chrom=as.factor(seg.seqnames),
                     maploc=seg.start, 
                     data.type="logratio",
                     sampleid=alt.bin.ids[ar.i['alt']]))
    smoothed.CNAo <- smooth.CNA(CNAo)
    seg.CNAo <- segment(smoothed.CNAo, verbose=1, alpha=0.01, eta=0.05, min.width=5)
    sd.D <- sd(D, na.rm=TRUE)
    
    number_of_chunks = ceiling(length(D) / 100)
    number_of_chunks=100
    sd.D <- sapply(split(seq_len(length(D)), 
                 cut(seq_len(length(D)), pretty(seq_len(length(D)), number_of_chunks))),
           function(x) sd(D[x], na.rm=TRUE))
    sd.D <- mean(sd.D, na.rm=TRUE)
    seg.CNAo$output$seg.sd <- sd.D
    
    seg.drift <- CCLid:::.estimateDrift(seg.CNAo, z.cutoff=1:3)
    seg.CNAo$output <- seg.drift$seg
    class(seg.CNAo) <- 'CCLid'
    CCLid:::plot.CCLid(seg.CNAo)
    return(list("frac"=seg.drift$frac,
                "cna.obj"=seg.CNAo))
  })
  names(cn.drift) <- alt.ref.idx$id
  save(cn.drift, file="~/cn_drift.rda")
  
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
