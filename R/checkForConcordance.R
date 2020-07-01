#' checkForConcordance
#' @description Checks for concordance between sample in first column of
#' input matrix and the rest of the reference matrix
#' 
#' @param x.mat Input matrix from compareVcf, or ref.mat
#' @param metric 'euclidean', 'jaccard', or 'cor'
#' @param meta.dat Metadata, defaults to meta.df if left NULL (default:NULL)
#' @param return.matrix Return distance matrix (Default: FALSE)
#' @param return.pred Return predicted Match/Nonmatch (Default: TRUE)
#' @param fill.na Fill NA's in the input matrix with median value (Default=TRUE)
#' @param rm.gcsi Remove GNE samples from the input matrix (Default=TRUE)
#' @param sampleID Sample ID (Default='Sample')
#' @param meta.df Cell line metadata, accessible from CCLid::ccl_table
#' @param verbose Verbose
#' @importFrom stats median
#' @importFrom utils data
#' 
#' @return A list containing Distance Matrix and list of Predictions (M/NM)
#' @export
#'
checkForConcordance <- function(x.mat, metric='euclidean', 
                                meta.dat=NULL, verbose=FALSE,
                                return.matrix=FALSE, return.pred=TRUE,
                                fill.na=TRUE, rm.gcsi=TRUE, sampleID='Sample',
                                meta.df){
  
  colnames(x.mat)[1] <- sampleID
  x.mat <- as.matrix(x.mat)
  storage.mode(x.mat) <- 'numeric'
  ret.dat <- list()
  
  if(fill.na | rm.gcsi){
    if(rm.gcsi){
      warning("Removing gCSI dataset from input matrix")
      gcsi.idx <- c(grep("GNE_", colnames(x.mat)),
                    grep("^UNK", colnames(x.mat)))
      x.mat <- x.mat[,-gcsi.idx]
    }
    if(fill.na){
      warning("Filling NA values in input matrix with median")
      
      x.mat[is.na(x.mat)] <- median(x.mat, na.rm =TRUE)
    }
  } else {
    warning("If input matrix has NA values, run-time is significantly longer")
  }
  
  x.dist <- CCLid::similarityMatrix(x.mat, method = metric)
  if(is.null(meta.dat)){
    if(verbose) print("Loading in metadata designed for reference matrix...")
    #data(meta.df)
    meta.dat <- meta.df
  }
  gc();
  
  if(return.pred){
    D.vals <- lapply(list("baf"=x.dist), CCLid::splitConcordanceVals, meta.df=meta.dat)
    balanced <- CCLid::balanceGrps(D.vals)
    models <- CCLid::trainLogit(balanced, predictors=c('baf'))
    pred <- assemblePredDat(D.vals, known.class=FALSE)
    gc();
    pred <- mkPredictions(pred, models)
    pred <- pred[order(pred$baf.fit),]
    pred$Var1 <- as.character(pred$Var1)
    pred$Var2 <- as.character(pred$Var2)
    gc();
    
    sample.idx <- grep(paste0("", sampleID, ""), pred$Var1, fixed = TRUE)
    sample.idx <- c(sample.idx, grep(paste0("", sampleID, ""), pred$Var2, fixed = TRUE))
    pred.sample <- pred[sort(sample.idx),]
    pred <- split(pred.sample, pred.sample$baf.p.fit)
    ret.dat[['pred']] <- pred
  }
  
  if(return.matrix) ret.dat[['D']] <- x.dist
  return(ret.dat)
}