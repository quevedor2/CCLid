.demoP <- function(data.type, n.sim=5, n.loci=1000, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  
  demo.dat <- switch(data.type,
         "BAF"={
           message("Generating randomized 'BAF' data...")
           mus <- c(0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9)
           sds <- rep(0.01, length(mus))
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
             p <- data.frame("0"=round(runif(n = n.loci, min=0, max=0.5),2),
                             "1"=round(runif(n = n.loci, min=0, max=0.5),2), 
                             check.names=F)
             p$`2` <- 1 - rowSums(p)
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

combineSamples <- function(){
  
}

demo <- genDemoData(data.type='BAF', n.sim=5, n.pop=10, n.loci=1000, seed=1)
test.sample <- demoSample('BAF', 1, demo$prob)
test.sample <- annoDemoMat(c(1), test.sample)

apply(demo$matrix, 2, function(i) cor(test.sample[,1], i))

corrplot::corrplot(cor(demo$matrix))
