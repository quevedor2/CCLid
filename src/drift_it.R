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
  alt.ds <- 'CCLE'
  
  vcf.dir <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs',
                       dataset)
  all.vcfs <- list.files(vcf.dir, pattern="vcf.gz$")
  rna.meta.df <- readinRnaFileMapping()
  names(all.vcfs) <- gsub(".snpOut.*", "", all.vcfs)
  
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
      sig.diff.gr <- lapply(x.drift$cna.obj, function(each.sample, sig.es=NULL){
        es <- each.sample$output
        sig.idx <- which(es$t >= 3)
        if(length(sig.idx) > 0){
          es <- es[sig.idx,] 
          sig.es <- lapply(split(es, es$ID), makeGRangesFromDataFrame, keep.extra.columns=TRUE)
        }
        return(sig.es)
      })
      frac.cnt <- x.drift$frac
      
      summ <- list("frac"=frac.cnt,
                   "sig"=sig.diff.gr)
    } else {
      summ=NULL
    }
    
    
    return(summ)
  }, mc.cores = 8)
  save(vcf.drift, file="~/vcf_drift.rda")

  .getDrift <- function(i, idx=1){
    if(length(i$frac) >= idx){
      delta <- i$frac[[idx]][3,,drop=FALSE]
      rownames(delta) <- gsub("RNA_", "", names(i$frac)[1])
    } else {
      delta <- NULL
    }

    if(!is.null(delta)) colnames(delta) <- gsub("_.*", "", colnames(delta))
    return(as.data.frame(delta))
  }

  ## Aggregate drift matrix
  rna.drift <- do.call(gtools::smartbind, lapply(vcf.drift, .getDrift, idx=1)) #RNA - GDSC/CCLE
  cl.drift <- do.call(gtools::smartbind, lapply(vcf.drift, .getDrift, idx=2)) # GDSC - CCLE
  colnames(rna.drift) <- paste0("RNA_", colnames(rna.drift))
  colnames(cl.drift) <- paste0("SNP_", colnames(cl.drift))
  rna.drift$ID <- rownames(rna.drift)
  cl.drift$ID <- rownames(cl.drift)
  rna.cl.drift <- merge(rna.drift, cl.drift, by="ID", all.y=TRUE)

  ## Measure overlap between RNA and cell lines
  seg.sig <- lapply(vcf.drift, function(i) if(length(i$sig) == 2) i$sig else NULL)
  null.idx <- which(sapply(seg.sig, is.null))

  # Interesct and fractionate
  all.drift.ovs <- lapply(seg.sig[-null.idx], function(seg){
    #seg <- seg.sig[-null.idx][[max.intersect]]
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

  # Assign the colours
  all.row.ids <- unique(as.character(sapply(idx.drift, rownames)))
  cols <- setNames(RColorBrewer::brewer.pal(length(all.row.ids), "Set2"),
                   all.row.ids)
  cols['intersect'] <- '#ffffcc'
  cols['no_drift'] <- 'grey'

  ## Plot amount of intra- and inter-institutional drift
  pdf("~/drift_datasets.pdf")
  rcl.sub.drift <- rna.cl.drift[match(names(all.drift.ovs), rna.cl.drift$ID),]
  par(mfrow=c(3,1), mar=c(0.5, 4.1, 0.5, 2.1))
  with(rcl.sub.drift[ord,], barplot(RNA_GDSC, ylab="RNA/GDSC", las=1, ylim=c(0,1), col=cols['RNA/GDSC'], xaxt='n', xlab=''))
  with(rcl.sub.drift[ord,], barplot(RNA_GDSC, ylab="GDSC/CCLE", las=1, ylim=c(0,1), col=cols['GDSC/CCLE'], xaxt='n', xlab=''))
  with(rcl.sub.drift[ord,], barplot(RNA_GDSC, ylab="RNA/CCLE", las=1, ylim=c(0,1), col=cols['RNA/CCLE'], xaxt='n', xlab=''))
  
  ## Plot concordance in drift estimates
  col.idx <- grep(paste0(alt.ds, "$"), colnames(rna.cl.drift))
  comp.drift <- rna.cl.drift[,col.idx,drop=FALSE]
  par(mfrow=c(3,1), mar=c(2, 4.1, 2, 2.1))
  plot(comp.drift, pch=16, col='black', xlim=c(0,1), ylim=c(0,1),
       xlab=paste0(dataset, " (", gsub("_.*", "", colnames(comp.drift)[1]), ") - ", alt.ds, " (SNP)"),
       ylab=paste0(dataset, " (", gsub("_.*", "", colnames(comp.drift)[2]), ") - ", alt.ds, " (SNP)"),
       main='Genetic drift between cell lines (Fraction of genome)')
  r <- round(cor(comp.drift, use="pairwise", method="pearson")[1,2], 2)
  text(c(1,1), labels = paste0("r = ", r, " (pearson)"), adj=1)
  abline(coef = c(0,1), col="grey", lty=2)

  ## Plot the drift overlap matrix
  par(mfrow=c(3,1), mar=c(2, 4.1, 2, 2.1))
  lapply(idx.drift, function(m){
    barplot(as.matrix(m[,ord]), xaxt='n', col=cols[rownames(m)])
    par(xpd=TRUE)
    legend(x = 1, y = 0, fill=cols[rownames(m)], 
           legend=rownames(m), horiz=TRUE, box.lwd = 0)
    par(xpd=FALSE)
  })
  
  ## Plot individual cell lines
  max.intersect <- names(which.max(idx.drift[[1]]['intersect',]))
  max.intersect2 <- colnames(idx.drift[[1]][,ord])[41]
  max.rnaRef <- rcl.sub.drift[which.max(rcl.sub.drift$RNA_GDSC),]$ID
  par(mfrow=c(3,1), mar=c(0.5, 4.1, 0.5, 2.1))
  multiDriftPlot(seg.sig[-null.idx][[max.intersect]], ref.ds=dataset, alt.ds=alt.ds)
  multiDriftPlot(seg.sig[-null.idx][[max.intersect2]], ref.ds=dataset, alt.ds=alt.ds)
  multiDriftPlot(seg.sig[-null.idx][[max.rnaRef]], ref.ds=dataset, alt.ds=alt.ds)
  dev.off()

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
