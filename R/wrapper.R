#' loadRef
#' @description Wrapper to load in all the reference matrices
#' and annotate them with groups, as well as calculate variance
#' @param pdir Path to directory containing or to download reference data
#' @param analysis Only support BAF at the moment
#' @param ... bin.size to group variance data
#' 
#' @return List containing "ref"=reference BAF matrix
#' and "var"=list containing variance for bin-sizes
#' @export
#'
#' @examples
#' PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/CCLid/CCLid"
#' loadRef(PDIR, 'baf', bin.size=5e5)
loadRef <- function(pdir=NULL, analysis='baf', ...){
  print("Attaching/downloading reference matrix...")
  ref.mat <- downloadRefCCL(toupper(analysis), saveDir = PDIR)
  
  print("Ensuring proper format and caluclating variance per bin...")
  format.dat <- formatRefMat(name=toupper(analysis), ref.mat=ref.mat, 
                             analysis=tolower(analysis), ...) #bin.size=5e5
  ref.mat <- format.dat$mat
  var.dat <- format.dat$var
  rm(format.dat)
  
  ## Assign group IDs (e.g. 22Rv1.cel -> GDSC_22Rv1)
  print("Assigning group IDs...")
  new.ids <- assignGrpIDs(ref.mat, meta.df)
  new.ids[duplicated(new.ids)] <- gsub("_", "2_",new.ids[duplicated(new.ids)])
  colnames(ref.mat) <- new.ids
  return(list("ref"=ref.mat,
              "var"=var.dat))
}

#' compareVcf
#' @description Checks a VCF file(s) against reference dataset
#' to look for similarity to any known cell lines
#' @param vcfFile Path to VCF file to check against reference datasets
#' @param var.dat Variance data (list)
#' @param ref.mat Reference matrix (matrix)
#'
#' @return Matrix: Containing reference matrix subsetted to common SNPs as 
#' the input VCF, as well as a left-joined VCF data
#' @export
#'
#' @examples
#'   ref.dat <- CCLid::loadRef(PDIR, 'baf', bin.size=5e5)
#'   ref.mat <- ref.dat$ref
#'   var.dat <- ref.dat$var
#'   vcfFile <- '/mnt/work1/users/home2/quever/xfer/A549.sample_id.vcf' ## A549 WES
#'   vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, ref.mat=ref.dat$ref)
compareVcf <- function(vcfFile, var.dat, ref.mat){
  vcf.map <- CCLid::mapVcf2Affy(vcfFile)
  
  ## Combine matrices and reduce features
  ## Find the overlap between the COMParator and the REFerence
  vcf.map.var <- CCLid::mapVariantFeat(vcf.map, var.dat)
  vcf.to.use <- vcf.map.var
  ov.idx <- CCLid::overlapPos(comp = vcf.to.use$BAF,
                              ref=ref.mat, mapping = 'probeset')
  x.mat <- cbind(vcf.to.use$BAF$BAF[ov.idx$comp], 
                 ref.mat[ov.idx$ref,])
  return(x.mat)
}


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
#' @param verbose 
#'
#' @return A list containing Distance Matrix and list of Predictions (M/NM)
#' @export
#'
#' @examples
#'   ref.dat <- CCLid::loadRef(PDIR, 'baf', bin.size=5e5)
#'   ref.mat <- ref.dat$ref
#'   var.dat <- ref.dat$var
#'   vcfFile <- '/mnt/work1/users/home2/quever/xfer/A549.sample_id.vcf' ## A549 WES
#'   vcf.mat <- compareVcf(vcfFile, var.dat=ref.dat$var, ref.mat=ref.dat$ref)
#'   checkForConcordance(vcf.mat)
checkForConcordance <- function(x.mat, metric='euclidean', 
                                meta.dat=NULL, verbose=FALSE,
                                return.matrix=FALSE, return.pred=TRUE,
                                fill.na=TRUE, rm.gcsi=TRUE, sampleID='Sample'){
  
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
    data(meta.df)
    meta.dat <- meta.df
  }
  
  if(return.pred){
    D.vals <- lapply(list("baf"=x.dist), CCLid::splitConcordanceVals, meta.df=meta.dat)
    balanced <- CCLid::balanceGrps(D.vals)
    models <- CCLid::trainLogit(balanced, predictors=c('baf'))
    pred <- assemblePredDat(D.vals, known.class=FALSE)
    pred <- mkPredictions(pred, models)
    pred <- pred[order(pred$baf.fit),]
    pred$Var1 <- as.character(pred$Var1)
    pred$Var2 <- as.character(pred$Var2)
    
    sample.idx <- grep(paste0("^", sampleID, "$"), pred$Var1)
    sample.idx <- c(sample.idx, grep(paste0("^", sampleID, "$"), pred$Var2))
    pred.sample <- pred[sort(sample.idx),]
    pred <- split(pred.sample, pred.sample$baf.p.fit)
    ret.dat[['pred']] <- pred
  }
  
  if(return.matrix) ret.dat[['D']] <- x.dist
  return(ret.dat)
}

