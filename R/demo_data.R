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
#' @param data.type Must be either "BAF" (default) or "geno" (case-sensitive), 
#' @param sample.ord Order
#' @param demo.dat Demo data
#'
#' @return
demoSample <- function(data.type='BAF', sample.ord, demo.dat){
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
#' $ meta  :'data.frame':	5 obs. of  5 variables:
#' ..$ ID: Factor w/ 5 levels "S1","S2","S3",..: 1 2 3 4 5
#' ..$ A : chr [1:5] "S1A" "S2A" "S3A" "S4A" ...
#' ..$ B : chr [1:5] "S1B" "S2B" "S3B" "S4B" ...
#' ..$ C : chr [1:5] "S1C" NA "S3C" "S4C" ...
#' ..$ D : chr [1:5] "S1D" NA NA NA ...
#' 
#' @param sample.ord Sample order
#' @param demo.mat Matrix of BAF/Geno
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
  meta.df <- reshape2::dcast(m.ids, ID ~ uniq)
  
  
  list("mat"=demo.mat, "meta"=meta.df)
}

#' genDemoData
#' @description Generates random demo data to play with
#'
#' @param data.type Must be "BAF" (default) or "geno" (case-sensitive)
#' @param n.pop Number of samples to simulate
#' @param ... 
#'
#' @return A list object
#'  matrix = A matrix of 'geno' or 'BAFs' for n.pop samples
#'  meta = a matrix of all by all sample comparisons
#'  prob = a probability matrix used to generate each sample
#'
#' @examples
#'  genDemoData(data.type='BAF', n.pop=5)
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
#' @param data.type 'BAF' or 'geno'
#' @param sample.mat Input matrix of BAF or geno values by samples
#' @param prop Proportion to mix the 2+ samples by
#'
#' @return
#' A weighted mean vector of BAF or a probability sampled genotype
#' @export
combineSamples <- function(data.type, sample.mat, prop){
  switch(data.type, 
         'BAF'=apply(sample.mat, 1, weighted.mean, w=prop),
         'geno'=apply(sample.mat, 1, sample, size=1, prob=prop))
}

.demo <- function(){
  data.type <- 'geno'
  data.type <- 'BAF'
  s.idx <- c(1,5)
  prop <- c(0.1, 0.9)
  
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

.demoDeconvolution <- function(){
  
  
  ## Test decomposition via NMF
  s.idx.nmf <- c(1,2,3)
  M <- demoSample(data.type, s.idx.nmf, demo$prob)
  M <- annoDemoMat(s.idx.nmf, M)$mat
  #M <- M[1:50,]
  
  prop.nmf <- c(0.9, 0.09, 0.01)
  prop.nmf <- c(0.60, 0.30, 0.10)
  A <- as.matrix(combineSamples('BAF', M, prop=prop.nmf))
  
  # Decomposition with complete data
  M1.mse <- .checkMse(A, M)
  M1 <-nnmf(A, k = 0, check.k = FALSE, init=list(W0 = M));
  as.matrix(M1$H[,1] / colSums(M1$H)) # [0.9, 0.09, 0.01]
  
  # Decomposition with complete data and a rank
  M1 <-nnmf(A, k = 1, check.k = FALSE, init=list(W0 = M));
  as.matrix(M1$H[,1] / colSums(M1$H)) # [0.003, 0.975, 0.021, 0]
  
  # Decomposition with incomplete data
  M2.mse <- .checkMse(A, M[,1:2])
  M2 <-nnmf(A, k = 0, check.k = FALSE, init=list(W0 = M[,1:2]));
  as.matrix(M2$H[,1] / colSums(M2$H)) # [0.905, 0.094]
  
  # Decomposition with incomplete data and rank
  M2 <-nnmf(A, k = 1, check.k = FALSE, init=list(W0 = M[,1:2]));
  as.matrix(M2$H[,1] / colSums(M2$H)) # [0.0001, 0.947, 0.05]
  
  # Decomposition with incomplete data
  M3.mse <- .checkMse(A, M[,2:3])
  M3 <-nnmf(A, k = 0, check.k = FALSE, init=list(W0 = M[,2:3]));
  as.matrix(M3$H[,1] / colSums(M3$H)) # [0.581, 0.418]
  
  # Decomposition with incomplete data and a rank
  M3 <-nnmf(A, k = 1, check.k = FALSE, init=list(W0 = M[,2:3]));
  as.matrix(M3$H[,1] / colSums(M3$H)) # [1, 0, 0]
  
  # Decomposition with incomplete data and noise 
  A2 <- A + runif(n=nrow(A), min = -0.2, max = 0.2)
  M4.mse <- .checkMse(A2, M)
  M4 <-nnmf(A2, k = 0, check.k = FALSE, init=list(W0 = M));
  as.matrix(M4$H[,1] / colSums(M4$H)) # []
  
}

.demoRna <- function(){

  ## Load in Ref mat file
  PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/CCLid/CCLid"
  ref.dat <- CCLid::loadRef(PDIR, 'baf', bin.size=5e5)
  
  ## Load in VCF file of external data
  vcfFile <- '/mnt/work1/users/home2/quever/xfer/A549.sample_id.vcf' ## A549 WES
  vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, ref.mat=ref.dat$ref)
  
  ## Look for similarity
  sample=basename(vcfFile)
  colnames(vcf.mat)[1] <- sample
  pred <- checkForConcordance(vcf.mat, sampleID=sample) 
  
  ## Look for drift
  all.ids <- unique(unlist(pred$pred$M[,c('Var1', 'Var2')]))
  bdf <- bafDrift(vcf.mat[,c(sample, all.ids[grep(sample, all.ids, invert = TRUE)])])
  
  ## Return finished object for WebApp
  list("pred"=pred$pred$M,
       "drift"=bdf$frac[[1]],
       "seg"=bdf$cna.obj[[sample]]$output)
}
