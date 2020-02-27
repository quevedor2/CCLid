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
  dataset <- 'GDSC'  # or GNE
  alt.ds <- 'CCLE'
  
  ## Get CIN scores
  drug.pset <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/PSets'
  psets <- loadInPSets(drug.pset)
  cin <- getCinScore(psets, 'sum')
  # merge
  cin.m <- as.data.frame(t(plyr::rbind.fill(lapply(cin, function(i) as.data.frame(t(i))))))
  colnames(cin.m) <- names(cin)
  cin.m$tCIN <- rowSums(cin.m)
  cin.m$pID <- rownames(cin.m)
  
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
  cn.z <- 1; b.z <- 5
  summ.frac <- summarizeFracDrift(cn.drifts=cn.drifts, cn.z=cn.z,
                                  baf.drifts=baf.drifts, baf.z=b.z,
                                  include.id=TRUE)
  drift.cin <- lapply(summ.frac, function(i){
    i$pID <- meta.df[match(i$ID, meta.df$ID),]$PharmacoGX_ID
    cin.i <- Reduce(function(x,y) merge(x,y,by='pID', all=TRUE), list("drift"=i, "cin"=cin.m, "abc"=abc))
    cin.i
  })
  cn.d <- drift.cin$cn
  baf.d <- drift.cin$baf
  
  cn.d.abc <- do.call(rbind, apply(cn.d[,-c(1:7)], 2, function(i){
    data.frame("cinR"=cor(cn.d$tCIN, i, use="complete.obs"),
               "driftR"=cor(cn.d$drift, i, use="complete.obs"))
  }))
  cn.d.abc[order(cn.d.abc$driftR),]
  
  baf.d.abc <- do.call(rbind, apply(baf.d[,-c(1:7)], 2, function(i){
    data.frame("cinR"=cor(baf.d$tCIN, i, use="complete.obs"),
               "driftR"=cor(baf.d$drift, i, use="complete.obs"))
  }))
  baf.d.abc[order(baf.d.abc$driftR),]
  
  
  
  cn.d <- cn.d[which(cn.d$drift > 0.1),]
  plot(cn.d$drift, cn.d$tCIN)
  cor(cn.d$drift, cn.d$tCIN, use = "complete.obs", method='spearman')
  plot(drift.cin$baf$drift, drift.cin$baf$tCIN)
  cor(drift.cin$baf$drift, drift.cin$baf$tCIN, use = "complete.obs")
  
}