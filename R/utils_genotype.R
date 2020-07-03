#' plotHist
#' @description Plots the histogram for matching/non-matching dataset
#'
#' @param D Returns from predictions, containins $M and $NM as elements of the list
#' @importFrom stats quantile
#' @importFrom graphics hist
#' @importFrom graphics lines
#' @importFrom stats density
#' @importFrom scales alpha
#' 
#' @export
plotHist <- function(D){
  nm.vals <- D$NM
  m.vals <- D$M
  
  all.vals <- c(nm.vals$value, m.vals$value)
  max.x <- ceiling(max(all.vals) + abs(quantile(all.vals, 0.05)))
  min.x <- floor(min(all.vals) - abs(quantile(all.vals, 0.05)))
  spc <- (max.x - min.x) / 40
  hist(m.vals$value, col=alpha("green",0.5), border='white', 
       breaks=seq(min.x, max.x, by=spc), prob = TRUE,  
       xlab = "euclidean-distance", main = "GDSC-CCLE")
  lines(density(m.vals$value), lwd = 2, col = "darkgreen")
  
  hist(nm.vals$value, col=alpha("grey", 0.50), border='white', 
       breaks=seq(min.x, max.x, by=spc), prob = TRUE, add=T)
  lines(density(nm.vals$value), lwd = 2, col = "black")
}

#### Private Functions ####
#' Melts dataframe
#' @importFrom reshape2 melt
#' @param m Square similarity matrix
#'
.meltDf <- function(m){
  diag(m) <- m[upper.tri(m)] <- NA
  melt.m <- melt(m)
  melt.m[-which(is.na(melt.m$value)),]
}

.createFormula <- function(predictors){
  lapply(predictors, function(p){
    as.formula(paste0("id ~ ", paste(p, collapse="+")))
  })
}

#' .genNonmatchSnpDiff
#' @description generates a matrix of euclidean distance between probesets of 
#' sample pairs if the sample id's dont' match (NON-MATCH)
#' @param col.len  Col length
#' @param dr.nm Nonmatch SNP
#' @importFrom matrixStats rowDiffs
.genNonmatchSnpDiff <- function(col.len, dr.nm){
  ## Initialize a matrix of random pairs
  r.idx <- matrix(sample(1:ncol(dr.nm), col.len*2, replace=TRUE), nrow=2)
  ids.match <- TRUE
  while(ids.match){
    ## Initialize a matrix of random pairs
    id.mat <- matrix(colnames(dr.nm)[r.idx], nrow=2)
    id.mat <- gsub("^.*?_", "", id.mat)
    
    m.idx <- apply(id.mat, 2, function(i) {i[1] == i[2]})
    ids.match <- any(m.idx) # check for matching pairs
    
    ## update matching pairs to new indices
    if(ids.match){
      r.idx[,which(m.idx)] <- sample(1:ncol(dr.nm), sum(m.idx)*2, replace=TRUE)
    }
  }
  
  d <- as.data.frame(apply(r.idx, 2, function(com){
    rowDiffs(as.matrix(dr.nm[,com]))
  }))
  rownames(d) <- rownames(dr.nm)
  colnames(d) <- apply(r.idx,2, function(j) paste(colnames(dr.nm)[j],collapse=':'))
  return(d)
}

#' z-statistic
#'
#' @param dat vector of numbers to calculate z stat from
#' @param p returns p-value instead of z-statistic
#' @importFrom stats sd
#' @importFrom stats pnorm
.zval <- function(dat, p=FALSE){
  z <- (dat - mean(dat)) / sd(dat)
  if(p) z <- 2*pnorm(-abs(z))
  return(z)
}