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

demoSample <- function(data.type, sample.ord, demo.dat){
  sapply(sample.ord, function(i){
    switch(data.type,
           "BAF"=rnorm(n=length(demo.dat[[i]]$mus), 
                       mean=demo.dat[[i]]$mus, 
                       sd=demo.dat[[i]]$sd),
           "geno"=apply(demo.dat[[i]], 1, sample, 
                        x = c(0:2), size = 1, replace=F))
    
  })
}

annoDemoMat <- function(sample.ord, demo.mat){
  s <- rep(1, length(sample.ord))
  while(any(duplicated(paste0(sample.ord, s)))){
    s <- s + duplicated(paste0(sample.ord, s))
  }
  colnames(demo.mat) <- paste0("S", paste0(sample.ord, toupper(letters[s])))
  rownames(demo.mat) <- paste0("SNP", c(1:nrow(demo.mat)))
  demo.mat
}

genDemoData <- function(data.type='BAF', n.pop=10, ...){
  ## Generate demo data matrix
  demo.dat <- .demoP(data.type=data.type, ...) ## probabilities
  sample.ord <- sample(c(1:length(demo.dat)), size=n.pop, replace=T)
  
  ## Create samples from probabilities
  demo.mat <- demoSample(data.type, sample.ord, demo.dat) 

  ## Assign column and row annotations
  demo.mat <- annoDemoMat(sample.ord, demo.mat)
  
  list("matrix"=demo.mat, "prob"=demo.dat)
}

combineSamples <- function(data.type, sample.mat, prop){
  switch(data.type, 
         'BAF'=apply(sample.mat, 1, weighted.mean, w=prop),
         'geno'=apply(sample.mat, 1, sample, size=1, prob=prop))
}

.demo <- function(){
  data.type <- 'geno'
  data.type <- 'BAF'
  s.idx <- c(1,5)
  prop <- c(0.2, 0.8)
  
  demo <- genDemoData(data.type=data.type, n.sim=20, 
                      n.pop=40, n.loci=1000, seed=1234, sd=0.1)
  
  test.sample <- demoSample(data.type, s.idx, demo$prob)
  test.sample <- annoDemoMat(s.idx, test.sample)
  if(length(s.idx) > 1){
    sample.x <- as.matrix(combineSamples(data.type, test.sample, prop=prop))
  } else {
    sample.x <- test.sample
  }
  
  concordance <- apply(demo$matrix, 2, function(i) {
    switch(data.type, 
           "BAF"=cor(sample.x[,1], i),
           "geno"=1-philentropy::distance(t(data.frame(sample.x[,1], i)), 
                                          method='jaccard'))
  })
  head(as.matrix(sort(concordance, decreasing = T)), n=15)
  
  sim <- apply(demo$matrix, 2, function(i){
    apply(demo$matrix, 2, function(j){
      switch(data.type,
             "BAF"=cor(i,j),
             "geno"=1-philentropy::distance(t(data.frame(i, j)), 
                                          method='jaccard'))
    })
  })
  sim[sim == 1] <- NA
  gplots::heatmap.2(sim, trace = 'none',  dendrogram = 'none', 
                    cexRow = 0.4, cexCol=0.4, key = FALSE)
}
