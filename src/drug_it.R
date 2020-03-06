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
  ## Get CIN scores
  drug.pset <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/PSets'
  psets <- loadInPSets(drug.pset)
  cin <- getCinScore(psets, 'mean')

  ### Select Dataset! ###
  cor.data<- list()
  dataset <- 'GNE'  # or GNE, GDSC
  alt.ds <- 'CCLE'
  
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
  if(dataset=='GDSC'){
    cn.z <- 1; b.z <- 5  # GDSC
  } else if (dataset=='GNE'){
    cn.z <- 1; b.z <- 2  # GNE
  }
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
  
  ## Check Variance of drift compared to CIN
  sapply(drift.cin, function(dat){
    bin <- seq(1,7,by=0.1)
    drift.by.cin <- split(dat$drift, cut(dat$tCIN, bin))
    df <- data.frame("var"=sapply(drift.by.cin, function(i) (max(i, na.rm=TRUE) - min(i, na.rm=TRUE))),
                     "cin"=bin[-1])
    df[df < 0] <- NA
    # with(df, plot(var, cin))
    with(df, cor(var, cin, use='complete.obs'))
  })
  #           cn       baf 
  # GDSC-CCLE 0.6590224 0.4862286 
  # GNE-CCLE  0.6714717 0.5304877
  
  
  ## Visualization
  # i <- cn.d[,grep("Afatinib", colnames(cn.d))]
  pdf(file.path(PDIR, "drug_it", 
                paste0(dataset, "-", alt.ds, "_cor-abc.pdf")), height=5, width=10)
  col.idx <- c(1:(ncol(cn.d) - ncol(abc) + 1))
  par(mfrow=c(1,3))
  plot(cn.d$drift, cn.d$tCIN, xlim=c(0, 1), ylim=c(0,8), col=scales::alpha("black", 0.5),
       xlab="Drift", ylab="mean CIN70", pch=16)
  points(baf.d$drift, baf.d$tCIN, col=scales::alpha("blue", 0.5), pch=15)
  legend("topright", pch=c(16,15), col=c("black", "blue"), box.lwd = 0,
         legend=c(paste0("CN  (r = ", round(cor(cn.d$drift, cn.d$tCIN, use = "complete.obs"),2), ")"), 
                  paste0("BAF (r = ", round(cor(baf.d$drift, baf.d$tCIN, use = "complete.obs"),2), ")")))
  cor.data[[dataset]][['cn']] <- corWithDrug(dat.d=cn.d, col.idx = col.idx, title='CN-drift', text.thresh=0.5)
  cor.data[[dataset]][['baf']] <- corWithDrug(dat.d=baf.d, col.idx = col.idx, title='BAF-drift', text.thresh=0.5)
  dev.off()
  scp.path <- "scp quever@192.168.198.99:"
  cat(paste0(scp.path, file.path(PDIR, "drug_it", paste0(dataset, "-", alt.ds, "_cor-abc.pdf")), ' .\n'))
  
  
  pdf(file.path(PDIR, "drug_it", 
                paste0("ALL_cor-abc.pdf")), height=6, width=6)
  ds1 <- cor.data$GDSC[[1]]
  ds2 <- cor.data$GNE[[1]]
  ds1$ID <- rownames(ds1); ds2$ID <- rownames(ds2)
  r.ds <- merge(ds1[,c('ID', 'cinR', 'driftR')],ds2[,c('ID', 'cinR', 'driftR')], by='ID')
  
  with(r.ds, plot(cinR.x, cinR.y, xlim=c(-0.3, 0.3), ylim=c(-0.3, 0.3), col=scales::alpha("black", 0.8),
                  xlab="Correlation (GDSC)", ylab="Correlation (GNE)", pch=18))
  l <- c(paste0("CIN:ABC  (r = ", round(cor(r.ds$cinR.x, r.ds$cinR.y, use = "complete.obs"),2), ")"))
  ds.col <- c('cn'='blue', 'baf'='red')
  for(d.type in c('cn', 'baf')){
    ds1 <- cor.data$GDSC[[d.type]]
    ds2 <- cor.data$GNE[[d.type]]
    ds1$ID <- rownames(ds1); ds2$ID <- rownames(ds2)
    r.ds <- merge(ds1[,c('ID', 'cinR', 'driftR')],ds2[,c('ID', 'cinR', 'driftR')], by='ID')
    
    with(r.ds, points(driftR.x, driftR.y, col=scales::alpha(ds.col[d.type], 0.7), pch=18))
    l <- c(l,  paste0(toupper(d.type), "-drift:ABC (r = ", round(cor(r.ds$driftR.x, r.ds$driftR.y, use = "complete.obs"),2), ")"))
  }
  legend("topright", pch=18, col=c("black", ds.col), box.lwd = 0, legend=l)
  dev.off()
  scp.path <- "scp quever@192.168.198.99:"
  cat(paste0(scp.path, file.path(PDIR, "drug_it", paste0("ALL_cor-abc.pdf")), ' .\n'))
  
  
  
  sapply(c("GDSC", "GNE"), function(ds){
    sapply(c("cn", "baf"), function(ty){
      cn <- cor.data[[ds]][[ty]]
      cn <- cn[order(cn$drift.q),]
      write.table(cn, file=file.path(PDIR, "drug_it", paste0(ds, "-", ty, "_table.csv")), quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
      cat(paste0(scp.path, file.path(PDIR, "drug_it", paste0(ds, "-", ty, "_table.csv")), ' .\n'))
      NA
    })
  })
  # gdsc <- rownames(cor.data[['GDSC']][[ty]])
  # gne <- rownames(cor.data[['GNE']][[ty]])
  
}