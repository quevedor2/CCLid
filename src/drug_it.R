###########################
#### Preliminary steps ####
###########################
## Run before executing any of the following 3 analyses
loadInData <- function(){
  library(PharmacoGx)
  library(CCLid)
  
  ## Set dataset colors
  dataset.cols <- setNames(RColorBrewer::brewer.pal(6, "Dark2"),
                           c("GDSC", "CCLE", "gCSI", "CGP", "Pfizer"))
  
  PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/CCLid/CCLid"
  ref.dat <- CCLid::loadRef(PDIR, 'baf', bin.size=5e5)
}

#############################################
#### Concordance between Drift and CIN70 ####
#############################################
## Compares BAF and CN drift scores with CIN70 estimates
## This part assumes that drift_it.R was run and the data
## was saved to disk to load in
iveCinItAlready <- function(){
  dataset <- 'GDSC'  # or GNE, GDSC
  alt.ds <- 'CCLE'
  
  ## Get CIN scores
  drug.pset <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/PSets'
  psets <- loadInPSets(drug.pset)
  cin <- getCinScore(psets, 'mean')
  # merge
  cin.m <- as.data.frame(t(plyr::rbind.fill(lapply(cin, function(i) as.data.frame(t(i))))))
  colnames(cin.m) <- names(cin)
  cin.m$tCIN <- rowMeans(cin.m[,c(dataset, alt.ds)], na.rm=TRUE)
  cin.m$pID <- rownames(cin.m)
  
  ## Get Gene Expr
  genes <- c('PTEN')
  gene.exprs <- getGeneExpr(psets, genes)
  exprs.m <- lapply(setNames(genes, genes), function(g){
    m <- as.data.frame(t(plyr::rbind.fill(lapply(gene.exprs, function(i) {
      idx <- match(g, rownames(i))
      as.data.frame(i[idx,,drop=FALSE])
    }))))
    colnames(m) <- names(gene.exprs)
    m$tCIN <- round(rowMeans(m[,c(dataset, alt.ds)], na.rm=TRUE),2)
    m$pID <- rownames(m)
    m[,c('tCIN', 'pID')]
  })
  
  
  
  ## Load drift data
  load(file=file.path(PDIR, "drift_it", 
                      paste0(dataset, "-", alt.ds, "_cn_drift.rda"))) #cn.drifts
  load(file=file.path(PDIR, "drift_it", 
                      paste0(dataset, "-", alt.ds, "_baf_drift.rda"))) #baf.drifts
  
  ## Load ABC drug sensitivity data
  abc.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/abc/rds'
  abcs <- readRDS(file.path(abc.dir, paste0(dataset, "-", alt.ds, ".rds")))
  abcs.m <- lapply(abcs, function(i) {
    abc <- as.data.frame((diag(i)))
    rownames(abc) <- colnames(i)
    abc
  })
  abc <- cbind(data.frame("pID"=colnames(abcs[[1]])), do.call(cbind, abcs.m))
  colnames(abc) <- c("pID", names(abcs))
  
  ## Compare drift to CIN70
  cn.z <- 1; b.z <- 5  # GDSC
  # cn.z <- 1; b.z <- 2  # GNE
  summ.frac <- summarizeFracDrift(cn.drifts=cn.drifts, cn.z=cn.z,
                                  baf.drifts=baf.drifts, baf.z=b.z,
                                  include.id=TRUE)
  # Append Gene Expr
  if(length(exprs.m) > 1) {
    exprs.mu <- Reduce(function(x,y) merge(x,y,by='pID', all=TRUE), exprs.m)
  } else {
    exprs.mu <- exprs.m[[1]][,c('pID', 'tCIN')]
  }
  colnames(exprs.mu) <- c('pID', genes)
  
  # Append CIN scores
  drift.cin <- lapply(summ.frac, function(i){
    i$pID <- meta.df[match(i$ID, meta.df$ID),]$PharmacoGX_ID
    cin.i <- Reduce(function(x,y) merge(x,y,by='pID', all=TRUE), list("drift"=i, "cin"=cin.m, "gene"=exprs.mu, "abc"=abc))
    cin.i
  })
  cn.d <- drift.cin$cn
  baf.d <- drift.cin$baf
  
  # i <- cn.d[,grep("Afatinib", colnames(cn.d))]
  col.idx <- c(1:(ncol(cn.d) - ncol(abc) + 1))
  corWithDrug(dat.d=cn.d, col.idx = col.idx)
  corWithDrug(dat.d=baf.d, col.idx = col.idx)

  cn.d.abc <- do.call(rbind, apply(cn.d[,-col.idx], 2, function(i){
    # plot(i, cn.d$tCIN)
    # plot(i, cn.d$drift)
    df <- data.frame("n"=rowSums(matrix(c(x,y), ncol=2)),
                     "cinR"=cor(cn.d$tCIN, i, use="complete.obs"),
                     # "cinR.p"=tryCatch({cor.test(cn.d$tCIN, i, use="complete.obs")$p.value}, error=function(e){NA}),
                     "driftR"=cor(cn.d$drift, i, use="complete.obs"),
                     "driftR.p"=tryCatch({cor.test(cn.d$drift, i, use="complete.obs")$p.value}, error=function(e){NA}))
    cbind(df, t(sapply(genes, function(g){cor(cn.d[,g], i, use='complete.obs')})))
  }))
  cn.d.abc$cin.q <- p.adjust(cn.d.abc$cinR.p, method="fdr")
  cn.d.abc$drift.q <- p.adjust(cn.d.abc$driftR.p, method="fdr")
  head(cn.d.abc[order(cn.d.abc$drift.q),],10)
  head(cn.d.abc[order(cn.d.abc$cin.q),],10)
  with(cn.d.abc, plot(cinR, driftR, pch=16, col=scales::alpha("black", 0.7)))
  abline(h=0, v=0)
  for(i in seq(0.1, 0.6, by=0.1)){
    with(cn.d.abc[which(cn.d.abc$cin.q  < i | cn.d.abc$drift.q < i),], points(cinR, driftR, col=scales::alpha("red", 0.3), pch=16))
  }
  
  
  baf.d.abc <- do.call(rbind, apply(baf.d[,-col.idx], 2, function(i){
    data.frame("cinR"=cor(baf.d$tCIN, i, use="complete.obs"),
               "cinR.p"=tryCatch({cor.test(cn.d$tCIN, i, use="complete.obs")$p.value}, error=function(e){NA}),
               "driftR"=cor(baf.d$drift, i, use="complete.obs"),
               "driftR.p"=tryCatch({cor.test(cn.d$drift, i, use="complete.obs")$p.value}, error=function(e){NA}))
  }))
  baf.d.abc$cin.q <- p.adjust(baf.d.abc$cinR.p, method="fdr")
  baf.d.abc$drift.q <- p.adjust(baf.d.abc$driftR.p, method="fdr")
  head(baf.d.abc[order(baf.d.abc$drift.q),],10)
  head(baf.d.abc[order(baf.d.abc$cin.q),],10)

  
  
  cn.d <- cn.d[which(cn.d$drift > 0.1),]
  plot(cn.d$drift, cn.d$tCIN)
  cor(cn.d$drift, cn.d$tCIN, use = "complete.obs", method='spearman')
  plot(drift.cin$baf$drift, drift.cin$baf$tCIN)
  cor(drift.cin$baf$drift, drift.cin$baf$tCIN, use = "complete.obs")
  
}