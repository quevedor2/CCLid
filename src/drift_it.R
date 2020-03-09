###########################
#### Preliminary steps ####
###########################
## Run before executing any of the following 3 analyses
loadInData <- function(){
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
  dataset <- 'GNE' #GDSC
  alt.ds <- 'CCLE' #CCLE
  
  ## Find variant features and isolate for cell lines shared in datasets
  ref.mat.var <- mapVariantFeat(ref.dat$ref, ref.dat$var)
  cl.idx <- CCLid::findCclPairs(meta.df, ref.mat.var, ds=c(alt.ds, dataset))
  m.cls.idx <- cl.idx[sapply(cl.idx, function(i) length(i) >= 2)]
  
  ## Get drift distance between CL pairs using BAF
  # ccl.id <- 'COR-L23';  x.mat=ref.mat.var; centering='median'
  # cl.pairs <- m.cls.idx[[ccl.id]]; ref.ds=dataset; segmenter='PCF'
  baf.drifts <- mclapply(m.cls.idx, CCLid::getBafDrifts, x.mat=ref.mat.var, 
                       ref.ds=dataset, alt.ds=alt.ds, segmenter='PCF', 
                       centering='none', mc.cores = 8)
  save(baf.drifts, file=file.path(PDIR, "drift_it", 
                             paste0(dataset, "-", alt.ds, "_baf_drift.rda")))
  load(file=file.path(PDIR, "drift_it", 
                      paste0(dataset, "-", alt.ds, "_baf_drift.rda")))
  
  ## Load in CN bins
  cn.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn/50kb_bins'
  cn.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn'
  cn.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnv_predictions/input/cn/log2r_bins'
  bins.file <- list.files(cn.dir, pattern="\\.bins\\.", include.dirs = FALSE)
  bins.file <- bins.file[grep(paste(c(dataset, alt.ds), collapse="|"), bins.file)]
  bins <- lapply(bins.file, function(b) readRDS(file.path(cn.dir, b)))
  names(bins) <- gsub("_.*", "",  bins.file)
  
  cn.drifts=CCLid::getCNDrifts(ref.l2r=assayData(bins[[dataset]]),
                               alt.l2r=assayData(bins[[alt.ds]]),
                               seg.id='exprs', raw.id='L2Rraw',
                               fdat=featureData(bins[[alt.ds]])@data,
                               cell.ids=names(baf.drifts), 
                               centering='extreme',
                               quantnorm=if(dataset=='GNE') TRUE else FALSE)
  # ref.l2r=assayData(bins[[dataset]]); alt.l2r=assayData(bins[[alt.ds]]);
  # seg.id='exprs'; raw.id='L2Rraw'; fdat=featureData(bins[[alt.ds]])@data;
  # cell.ids=names(baf.drifts); segmenter='PCF';centering='extreme'
  save(cn.drifts, file=file.path(PDIR, "drift_it", 
                                paste0(dataset, "-", alt.ds, "_cn_drift.rda")))
  

  ## Load in data
  load(file=file.path(PDIR, "drift_it", 
                      paste0(dataset, "-", alt.ds, "_cn_drift.rda")))
  load(file=file.path(PDIR, "drift_it", 
                      paste0(dataset, "-", alt.ds, "_baf_drift.rda")))
  
  if(dataset == 'GDSC'){ cn.z <- 1; b.z <- 5 } else if(dataset == 'GNE'){  cn.z <- 1; b.z <- 2  }
  summ.frac <- summarizeFracDrift(cn.drifts=cn.drifts, cn.z=cn.z,
                                  baf.drifts=baf.drifts, baf.z=b.z,
                                  include.id=TRUE)
  pdf(file=file.path(PDIR, "drift_it", paste0(dataset, "-", alt.ds, "_baf-cn-frac.pdf")),
      width=5, height=5)
  plotFracDrift(summ.frac)
  dev.off()
  cat(paste0("scp quever@192.168.198.99:", 
             file.path(PDIR, "drift_it", paste0(dataset, "-", alt.ds, "_baf-cn-frac.pdf .\n"))))
  
  
  
  ## Find overlap between GRanges of BAF to CN drift
  df.cn <- cn.drifts$cna.obj$output
  gr.cn <- makeGRangesFromDataFrame(df.cn, keep.extra.columns = TRUE)
  gr.cn <- split(gr.cn, gr.cn$ID)
  df.baf <- lapply(baf.drifts, function(i) {
    d <- i$sig.gr[[1]]$output
    d$ID <- gsub("(GNE_)|(GDSC_)|(CCLE_)", "", names(i$sig.gr)[1])
    return(d)
  })
  gr.baf <- makeGRangesFromDataFrame(do.call(rbind, df.baf), keep.extra.columns = TRUE)
  gr.baf <- split(gr.baf, gr.baf$ID)
  
  
  ## Calculate concordance between SNP and CN drifted regions
  baf.drift.dat <- mclapply(c(0:9), function(baf.z){
    print(baf.z)
    driftOverlapMetric(gr.baf = gr.baf, gr.cn = gr.cn, 
                       cell.ids = names(baf.drifts),
                       baf.z=baf.z, cn.z=cn.z, cn.gtruth=FALSE)
  }, mc.cores = 10)
  drift.dat <-  driftOverlapMetric(gr.baf = gr.baf, gr.cn = gr.cn, 
                                   cell.ids = names(baf.drifts),
                                   baf.z=b.z, cn.z=cn.z, cn.gtruth=FALSE)
  
  ## Plot the saturation-sensitvity curve
  pdf(file=file.path(PDIR, "drift_it", paste0(dataset, "-", alt.ds, "_baf-cn-drift.pdf")),
      width=4, height=4)
  lapply(baf.drift.dat, function(drift.dat){
    with(drift.dat$sens, plot(x=x, y=y, pch=16, las=1, col=scales::alpha("grey", 0.8),
                              main="Agreement between inference of CN and BAF drift", xaxt='n',
                              ylim=c(0,1), ylab="Sensitivity", xlab="Conc. Threshold"))
    axis(side = 1, at = seq(0,1, by=0.2), labels = rev(seq(0,1, by=0.2)), las=1)
    
    ## add saturation points
    sat.points <- as.character(c(0.7, 0.8, 0.9))
    sat.cols <- c("#feb24c", "#fd8d3c", "#f03b20")
    
    tryCatch({
      lines(drift.dat$sens$x, predict(drift.dat$model),lty=2,col="black",lwd=2)
      abline(v=drift.dat$saturation[sat.points,], col=sat.cols)
      text(x =drift.dat$saturation[sat.points,], y=rep(1, length(sat.points)), 
           labels = sat.points, cex=0.8, adj=0 )
      print(1-drift.dat$saturation[sat.points,])
    }, error=function(e){ NA })
    #GNE(b.z=2) 0.7   0.8   0.9 
    #GNE(b.z=2) 0.867 0.823 0.747 
  })
  dev.off()
  cat(paste0("scp quever@192.168.198.99:", file.path(PDIR, "drift_it", paste0(dataset, "-", alt.ds, "_baf-cn-drift.pdf .\n"))))
  
  # baf.idx <- 4
  # as.matrix(head(sort(colSums(baf.drift.dat[[baf.idx]]$dat)), 30))
  # tail(sort(colSums(baf.drift.dat[[baf.idx]]$dat)), 30)
  
  ## Plot the drift CN-BAF examples
  # ccl.id <-'IGR-37'  #'786-0', 'HT-29', 'CL-40', 'SW1463', 'A172', 'MCAS', 'HCC1937', 'Namalwa', 'PC-3'
  # sapply(names(head(sort(colSums(drift.dat$dat)), 10)), function(ccl.id){
#   sapply(c('KE-37',
#            'COR-L23', 'JHOS-2', 'HCC1937', 'RS4-11',
# #           'NCI-H2029', 'HLE', 'HCC-366',
# #           'HCC1937', 'HuP-T4', 'COR-L23',
# #           'KNS-62', 'NCI-H522', 'CAS-1',
# #           "SW403", "VM-CUB-1", "HuH-6"
# #           'JHOS-2', 'NB-1', 'NCI-H23', 'PC-3', 
# #           'AsPC-1', 'A172', 
#            'Capan-1'), function(ccl.id){
  #'NCI-H23', 'NCI-H2052', 'MCF-7', 'A3-KAW'
  sapply(c('HCC1599', 'Hs-294T'), function(ccl.id){
    pdf(file=file.path(PDIR, "drift_it", paste0(dataset, "-", alt.ds, "_baf-cn-drift_", ccl.id, ".pdf")),
        width=7, height=3)
    CNAo <- cn.drifts$cna.obj
    CNAo$output <- split(CNAo$output, CNAo$output$ID)[[ccl.id]]
    CNAo$data <- CNAo$data[,c(1,2,grep(paste0("^", ccl.id, "$"), colnames(CNAo$data)))]
    CCLid:::plot.CCLid(CNAo, min.z=1)
    CCLid:::plot.CCLid(baf.drifts[[ccl.id]]$sig.gr[[1]], min.z=b.z)
    dev.off()
    
    meta.cclid <- meta.df[grep(paste0("^", ccl.id, "$"), meta.df$ID),]
    scp.path <- "scp quever@192.168.198.99:"
    path.tmp <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
    cat(paste0(scp.path, file.path(path.tmp, "GNE", "eacon", meta.cclid$GNE, "ASCAT", "L2R", "*png "), paste0("GNE_", meta.cclid$ID, ".png\n")))
    cat(paste0(scp.path, file.path(path.tmp, "CCLE", "eacon", meta.cclid$CCLE, "ASCAT", "L2R", "*png "), paste0("CCLE_", meta.cclid$ID, ".png\n")))
    cat(paste0(scp.path, file.path(path.tmp, "GDSC", "eacon", gsub(".cel", "", meta.cclid$GDSC, ignore.case=TRUE), "ASCAT", "L2R", "*png "), paste0("GDSC_", meta.cclid$ID, ".png\n")))
    cat(paste0(scp.path, file.path(PDIR, "drift_it", paste0(dataset, "-", alt.ds, "_baf-cn-drift_", ccl.id, ".pdf .\n"))))
  })
}

###################################
#### Drift between SNP and RNA ####
###################################
## This process is meant to genetic fraction
## of drift found in the RNAseq and the SNP array
## data.  As well as calculate the overlap of
## segments between the two
driftTech <- function(){
  # dataset <- 'GDSC'
  # alt.ds <- 'CCLE'
  dataset <- 'CCLE'
  alt.ds <- 'GDSC'
  
  vcf.dir <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs',
                       dataset)
  all.vcfs <- list.files(vcf.dir, pattern="vcf.gz$")
  rna.meta.df <- readinRnaFileMapping()
  names(all.vcfs) <-  sapply(gsub(".snpOut.*", "", all.vcfs), function(i){
    idx <- switch(dataset,
                  "GDSC"=grep(paste0("^", i, "$"), rna.meta.df$EGAF),
                  "CCLE"=grep(paste0("^", i, "$"), rna.meta.df$SRR))
    if(length(idx) >= 1){
      rna.meta.df[idx,]$ID[1]
    } else {
      gsub(".snpOut.vcf.gz$", "", i)
    }
  })
  vcf.ids <- setNames(names(all.vcfs), all.vcfs)
  
  ## Compare every VCF to the entire ref matrix to calculate BAF drift
  vcf.drift <- list()
  # ccl.id <- 'HLE'
  # vcf <- all.vcfs[ccl.id]
  # vcfFile=file.path(vcf.dir, vcf)
  for(vcf in all.vcfs){
    vcf.drift[[vcf]] <- getVcfDrifts(vcfFile=file.path(vcf.dir, vcf), 
                                     ref.dat, rna.meta.df, min.depth=5)
  }
  ## Very weird bugs happens when I use lapply
  # vcf.drift <- mclapply(all.vcfs[1:4], function(vcf){  
  #   getVcfDrifts(vcfFile=file.path(vcf.dir, vcf), ref.dat, rna.meta.df)
  # }, mc.cores = 3)
  names(vcf.drift) <- vcf.ids[names(vcf.drift)]
  save(vcf.drift, file=file.path(PDIR, "drift_it", 
                                 paste0(dataset, "-", alt.ds, "_vcf_drift.rda")))

  
  load(file=file.path(PDIR, "drift_it", 
                      paste0(dataset, "-", alt.ds, "_vcf_drift.rda")))
  load(file=file.path(PDIR, "drift_it", 
                      paste0(dataset, "-", alt.ds, "_baf_drift.rda")))
  
  ##################################################################################################
  gr.rna <- lapply(setNames(c("GDSC", "CCLE"),c("GDSC", "CCLE")), function(ds){
    rna <- lapply(vcf.drift, function(i) {
      d <- i$cna.obj[[1]]$output
      d <- d[grep(paste0(ds, "_"), d$ID),]
      d$ID <- gsub("(GDSC_)|(CCLE_)", "", d$ID)
      return(d)
    })
    rna <- rna[-which(sapply(rna, nrow)==0)]
    grna <- makeGRangesFromDataFrame(do.call(rbind, rna), keep.extra.columns = TRUE)
    grna <- split(grna, grna$ID)
    return(grna)
  })
  
  df.baf <- lapply(baf.drifts, function(i) {
    d <- i$sig.gr[[1]]$output
    d$ID <- gsub("(GDSC_)|(CCLE_)", "", names(i$sig.gr)[1])
    return(d)
  })
  gr.baf <- makeGRangesFromDataFrame(do.call(rbind, df.baf), keep.extra.columns = TRUE)
  gr.baf <- split(gr.baf, gr.baf$ID)
  
  
  
  
  rna.drift.dat <- mclapply(c(1:5), function(baf.z){
    print(baf.z)
    driftOverlapMetric(gr.baf = gr.baf, gr.cn = gr.rna$CCLE, 
                       cell.ids = names(baf.drifts),
                       baf.z=5, cn.z=baf.z)
  }, mc.cores = 5)

  pdf(file=file.path(PDIR, "drift_it", paste0(dataset, "-", alt.ds, "_baf-rna-drift.pdf")),
      width=5, height=5)
  lapply(rna.drift.dat, function(drift.dat){
    with(drift.dat$sens, plot(x=x, y=y, pch=16, las=1, col=scales::alpha("grey", 0.8),
                              main="Agreement between inference of CN and BAF drift", xaxt='n',
                              ylim=c(0,1), ylab="Sensitivity", xlab="Conc. Threshold"))
    axis(side = 1, at = seq(0,1, by=0.2), labels = rev(seq(0,1, by=0.2)), las=1)
    lines(drift.dat$sens$x, predict(drift.dat$model),lty=2,col="black",lwd=2)
    ## add saturation points
    sat.points <- as.character(c(0.7, 0.8, 0.9))
    sat.cols <- c("#feb24c", "#fd8d3c", "#f03b20")
    abline(v=drift.dat$saturation[sat.points,], col=sat.cols)
    text(x =drift.dat$saturation[sat.points,], y=rep(1, length(sat.points)), 
         labels = sat.points, cex=0.8, adj=0 )
    print(1-drift.dat$saturation[sat.points,])
    # 0.7     0.8     0.9 
    # 0.838   0.784   0.690 
  })
  dev.off()
  cat(paste0("scp quever@192.168.198.99:", file.path(PDIR, "drift_it", paste0(dataset, "-", alt.ds, "_baf-rna-drift.pdf .\n"))))
  
  head(sort(colSums(rna.drift.dat[[3]]$dat)), 50)
  tail(sort(colSums(rna.drift.dat[[3]]$dat)), 50)
  
  pdf("~/test2.pdf")
  sapply(c('VM-CUB-1', 'KM-H2', 'CL-40', 'HLE'), function(ccl.id){
    print(length(baf.drifts[[ccl.id]]))
    print(length(vcf.drift[[ccl.id]]))
    CCLid:::plot.CCLid(baf.drifts[[ccl.id]]$sig.gr[[1]], min.z=5)
    CCLid:::plot.CCLid(vcf.drift[[ccl.id]]$cna.obj[[1]], min.z=3)
  })
  dev.off()
  
  
  ## Drift overlaps between a given sample
  rna.meta.df[match(names(seg.sig[-null.idx])[cl.idx], rna.meta.df$EGAF),]
  ##################################################################################################
  ccl.id <- "CL-40"
  srr.file <- as.character(rna.meta.df[grep(ccl.id, rna.meta.df$ID),]$SRR)
  egaf.file <- as.character(rna.meta.df[grep(ccl.id, rna.meta.df$ID),]$EGAF)
  # col.idx <- grep("(UNK)|(GNE_)", colnames(ref.dat$ref))
  # ref.dat$ref <- ref.dat$ref[,-col.idx]
  vcf.x <- getVcfDrifts(vcfFile = file.path(vcf.dir, paste0(egaf.file, ".snpOut.vcf.gz")), 
                        ref.dat=ref.dat, rna.meta.df = rna.meta.df)
  multiDriftPlot(vcf.x$sig, chr.size.gr=CCLid:::.getChrLength(), ref.ds=dataset, alt.ds=alt.ds)
  
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
    multiDriftPlot(seg.sig[-null.idx][[cl.ids]], chr.size.gr=.getChrLength(), 
                   ref.ds=dataset, alt.ds=alt.ds)
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