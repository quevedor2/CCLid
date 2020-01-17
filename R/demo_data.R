.demoP <- function(data.type, n.sim=5, n.loci=1000, seed=NULL, sd=0.01){
  if(!is.null(seed)) set.seed(seed)
  
  demo.dat <- switch(data.type,
         "BAF"={
           message("Generating randomized 'BAF' data...")
           mus <- c(0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9)
           sds <- rep(sd, length(mus))
           lapply(c(1:n.sim), function(i){
             components <- sample(c(1:length(mus)), 
                                  prob=c(0.25, 0.025, 0.05, 0.35, 0.05, 0.025, 0.25),
                                  size=n.loci, replace=TRUE) 
             list("mus"=mus[components], "sd"=sds[components])
           })
         },
         "geno"=  {
           message("Generating randomized 'Genotype' data...")
           lapply(c(1:n.sim), function(i) {
             r <- rbinom(n = n.loci, size = 2, prob = 0.5)
             p <- matrix(rep(0.10, n.loci*3), ncol=3)
             for(i in seq_along(r)){
               p[i,r[i]] <- 0.8
             }
             p
           })
         })
}

#' demoSample
#' @description Randomize a sample BAF or Geno giving probabilities
#' 
#' @param data.type 
#' @param sample.ord 
#' @param demo.dat 
#'
#' @return
#' @export
#'
#' @examples
demoSample <- function(data.type, sample.ord, demo.dat){
  sapply(sample.ord, function(i){
    switch(data.type,
           "BAF"=rnorm(n=length(demo.dat[[i]]$mus), 
                       mean=demo.dat[[i]]$mus, 
                       sd=demo.dat[[i]]$sd),
           "geno"=apply(demo.dat[[i]], 1, function(j){
              sample(x=c(0:2), size=1, prob=j, replace=FALSE)
             }))
    
  })
}

#' annoDemoMat
#' @description Annotate the demo SNP x Sample matrix to generate meta.df
#' 
#' @param sample.ord 
#' @param demo.mat 
#'
#' @return
#' @export
annoDemoMat <- function(sample.ord, demo.mat){
  s <- rep(1, length(sample.ord))
  while(any(duplicated(paste0(sample.ord, s)))){
    s <- s + duplicated(paste0(sample.ord, s))
  }
  colnames(demo.mat) <- paste0("S", paste0(sample.ord, toupper(letters[s])))
  rownames(demo.mat) <- paste0("SNP", c(1:nrow(demo.mat)))
  
  ## Create the meta.df data structure
  m.ids <- data.frame("ID"=gsub("[A-Z]$", "", colnames(demo.mat)),
                      "uniq"=gsub("^S[0-9]*", "", colnames(demo.mat)),
                      "IDs"=colnames(demo.mat))
  meta.df <- dcast(m.ids, ID ~ uniq)
  
  
  list("mat"=demo.mat, "meta"=meta.df)
}

#' genDemoData
#' @description Generates random demo data to play with
#'
#' @param data.type 
#' @param n.pop 
#' @param ... 
#'
#' @return
#' @export
genDemoData <- function(data.type='BAF', n.pop=10, ...){
  ## Generate demo data matrix
  demo.dat <- .demoP(data.type=data.type, ...) ## probabilities
  sample.ord <- sample(c(1:length(demo.dat)), size=n.pop, replace=T)
  
  ## Create samples from probabilities
  demo.mat <- demoSample(data.type, sample.ord, demo.dat) 

  ## Assign column and row annotations
  demo <- annoDemoMat(sample.ord, demo.mat)
  
  list("matrix"=demo$mat, "meta"=demo$meta, "prob"=demo.dat)
}

#' combineSamples
#' @description Combines the BAF or Geno or two or samples given a set proportion
#' 
#' @param data.type 
#' @param sample.mat 
#' @param prop 
#'
#' @return
#' @export
combineSamples <- function(data.type, sample.mat, prop){
  switch(data.type, 
         'BAF'=apply(sample.mat, 1, weighted.mean, w=prop),
         'geno'=apply(sample.mat, 1, sample, size=1, prob=prop))
}

.demo <- function(){
  library(CCLid)
  data.type <- 'geno'
  data.type <- 'BAF'
  s.idx <- c(1,5)
  prop <- c(0.5, 0.5)
  
  demo <- genDemoData(data.type=data.type, n.sim=20, 
                      n.pop=40, n.loci=1000, seed=1234, sd=0.1)
  
  meta.df <- demo$meta
  test.sample <- demoSample(data.type, s.idx, demo$prob)
  test.sample <- annoDemoMat(s.idx, test.sample)$mat
  if(length(s.idx) > 1){
    sample.x <- as.matrix(combineSamples(data.type, test.sample, prop=prop))
  } else {
    sample.x <- test.sample
  }
  
  

  d.mat <- similarityMatrix(demo$matrix, 'euclidean')
  new.ids <- assignGrpIDs(d.mat, meta.df)
  colnames(d.mat) <- rownames(d.mat) <- as.character(new.ids)
  
  D.vals <- lapply(list("baf"=d.mat), splitConcordanceVals, meta.df=meta.df)
  balanced <- balanceGrps(D.vals)
  plotHist(D.vals$baf)
  
  models <- trainLogit(balanced, predictors=c('baf'))
  
  pred <- assemblePredDat(D.vals, known.class=TRUE)
  pred <- mkPredictions(pred, models)
  
  ## Viz
  head(pred[order(pred$baf.fit),])
  ggplot(pred,aes(x=baf,y=baf.fit)) + stat_binhex()
  
  
  
  ## Compare sampleX to ref data
  x.mat <- cbind(sample.x, demo$matrix)
  x.dist <- similarityMatrix(x.mat, 'euclidean')[,1,drop=FALSE]

  x.vals <- lapply(list("baf"=x.dist), splitConcordanceVals, meta.df=NULL)
  pred <- assemblePredDat(x.vals, known.class=FALSE)
  pred <- mkPredictions(pred, models)
  head(pred[order(pred$baf.fit),], n=15)
  ggplot(pred,aes(x=baf,y=baf.fit)) + stat_binhex()
  
  
  sim <- d.mat 
  sim[sim == 1] <- NA
  gplots::heatmap.2(sim, trace = 'none',  dendrogram = 'none', 
                    cexRow = 0.4, cexCol=0.4, key = FALSE)
}

.demoRna <- function(){
  library(VariantAnnotation)
  library(CCLid)
  
  PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/CCLid"
  analysis <- 'baf'
  ref.mat <- downloadRefCCL("BAF", saveDir = PDIR)
  ref.mat <- formatRefMat(ref.mat, analysis='baf')
  
  vcfFile='/mnt/work1/users/home2/quever/xfer/NGSST_04.mutect2.hg19.vcf' ## Mix
  vcfFile='/mnt/work1/users/home2/quever/xfer/NGSST_05.mutect2.hg19.vcf' ## Mix
  vcfFile <- '/mnt/work1/users/home2/quever/xfer/A549.sample_id.vcf' ## A549 WES
  vcf.map <- mapVcf2Affy(vcfFile)
  
  ## Find the overlap between the COMParator and the REFerence
  ov.idx <- overlapPos(comp = vcf.map$BAF,
                       ref=ref.mat, mapping = 'probeset')
  
  ## Combine matrices and reduce features
  x.mat <- cbind(vcf.map$BAF$BAF[ov.idx$comp], 
                 ref.mat[ov.idx$ref,])
  ref.vars <- apply(ref.mat, 1, var)
  
  
  ## Calculate evidence of genetic drift between sample X and A549
  ccle.id <- 'ARLES_p_NCLE_DNAAffy2_S_GenomeWideSNP_6_B06_256036'
  gdsc.id <- '^A549$'
  x.drift <- bafDrift(sample.mat=x.mat[,c(grep(ccle.id, colnames(x.mat)), 
                                          grep(gdsc.id, colnames(x.mat)), 
                                          1, 5)])
  plot(x.drift$cna.obj[[ccle.id]])
  
  
  ## Look for overall similarity

  # x.mat[is.na(x.mat)] <- median(as.matrix(x.mat), na.rm=T)
  x.dist <- similarityMatrix(x.mat, 'cor')[,1,drop=FALSE]
  as.matrix(c(head(x.dist[order(x.dist, decreasing = TRUE),], 10),
              head(x.dist[order(x.dist, decreasing = FALSE),], 10)))
  
  ccle.id <- 'ARLES_p_NCLE_DNAAffy2_S_GenomeWideSNP_6_B06_256036'
  gdsc.id <- '^A549$'
  gdsc.id <- '^DU-145$'
  gdsc.id <- '^HCC2998$'
  
  x.mat2 <- x.mat[,c(grep(ccle.id, colnames(x.mat)), grep(gdsc.id, colnames(x.mat)))]
  x.mat2 <- x.mat[,c(1, grep(gdsc.id, colnames(x.mat)))]
  plot(x.mat2, xlim=c(-0.5, 1.5), ylim=c(-0.5, 1.5))
}