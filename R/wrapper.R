#' loadRef
#' @description Wrapper to load in all the reference matrices
#' and annotate them with groups, as well as calculate variance
#'
#' @param analysis Only support BAF at the moment
#' @param ... Extra param
#' @param PDIR Directory for download data
#' @param rm.gne Remove gCSI from the analysis (default=FALSE)
#' @param bin.size Bin size (default=1e6)
#' @param verbose Verbose
#' 
#' @return List containing "ref"=reference BAF matrix
#' and "var"=list containing variance for bin-sizes
#' @export
#'
loadRef <- function(PDIR=NULL, analysis='baf', rm.gne=FALSE, bin.size=1e6, verbose=FALSE, ...){
  # if(verbose) print("Checking for existing reference data...")
  # ref.data.exists <- any(grepl(paste0(as.integer(bin.size), ".", toupper(analysis)), list.files(PDIR)))
  # 
  # if(ref.data.exists){
  #   if(verbose) print(paste0("Reading existing data: ", paste0(as.integer(bin.size), ".", toupper(analysis), ".rds"), "..."))
  #   ref.dat <- file.path(PDIR, paste0(as.integer(bin.size), ".", toupper(analysis), ".rds"))
  #   format.dat <- readRDS(ref.dat)
  # } else {
  # 
  # }
  if(verbose) print("Attaching/downloading reference matrix...")
  ref.mat <- downloadRefCCL(toupper(analysis), saveDir = PDIR, verbose=verbose)
  
  if(rm.gne){
    rm.idx <- grep("^Unk*", colnames(ref.mat))
    ref.mat <- ref.mat[,-rm.idx]
  }
  
  if(verbose) print("Ensuring proper format and caluclating variance per bin...")
  format.dat <- formatRefMat(name=toupper(analysis), ref.mat=ref.mat, saveDir = PDIR, 
                             analysis=tolower(analysis), bin.size=bin.size, ...) #bin.size=5e5
  
  ref.mat <- format.dat$mat
  var.dat <- format.dat$var
  rm(format.dat)
  
  ## Assign group IDs (e.g. 22Rv1.cel -> GDSC_22Rv1)
  if(verbose) print("Assigning group IDs...")
  if(!file.exists(file=file.path(PDIR, "col_ids.rda"))){
    new.ids <- assignGrpIDs(ref.mat, meta.df)
    new.ids[duplicated(new.ids)] <- gsub("_", "2_",new.ids[duplicated(new.ids)])
    save(new.ids, file=file.path(PDIR, "col_ids.rda"))
  } else {
    load(file=file.path(PDIR, "col_ids.rda"))
  }
  colnames(ref.mat) <- new.ids
  return(list("ref"=ref.mat,
              "var"=var.dat))
}

#' compareVcf
#' @description Checks a VCF file(s) against reference dataset
#' to look for similarity to any known cell lines
#'
#' @param vcfFile Path to VCF file to check against reference datasets
#' @param var.dat Variance data (list)
#' @param ref.mat Reference matrix (matrix)
#' @param max.snps Max number of SNPs to reduce 
#' @param ids IDs
#' @param sampletype Strictly for labelling purposes
#' @param ... Extra param
#' @param all.ids IDs to subset to
#'
#' @return Matrix: Containing reference matrix subsetted to common SNPs as 
#' the input VCF, as well as a left-joined VCF data
#' @export
#'
compareVcf <- function(vcfFile, var.dat, ref.mat, 
                       max.snps=1e6, ids=NULL, sampletype='RNA', all.ids=NULL, ...){
  vcf.map <- CCLid::mapVcf2Affy(vcfFile)
  vcf.map <- .filt(vcf.map, ...) ## Memory: up to 1.8Gb 
  
  ## Combine matrices and reduce features
  ## Find the overlap between the COMParator and the REFerence
  vcf.map.var <- CCLid::mapVariantFeat(vcf.map, var.dat)
  vcf.to.use <- vcf.map.var
  ov.idx <- CCLid::overlapPos(comp = vcf.to.use$BAF,
                              ref=ref.mat, mapping = 'probeset')
  if(nrow(ov.idx) > max.snps){
    ov.idx <- ov.idx[order(ov.idx$ref)[1:max.snps],]
  }
  rm(vcf.map, vcf.map.var); gc()  ## Cleanup
  
  ## BIG memory sink: shoots up to 7Gb
  if(is.null(ids)){
    refm <- ref.mat
  } else {
    colidx <- which(colnames(ref.mat) %in% all.ids)
    message("Isolating for: ", paste(colnames(ref.mat)[colidx], collapse=", "))
    refm <- ref.mat[,colidx,drop=FALSE]
  }
  x.mat <- cbind(vcf.to.use$BAF$BAF[ov.idx$comp], 
                 refm[ov.idx$ref,])
  colnames(x.mat)[1] <- paste0(sampletype, "_", gsub(".vcf.*", "", basename(vcfFile)))
  
  if(storage.mode(ref.mat[,1]) == 'integer'){
    x.mat[,-1] <- x.mat[,-1] / 100
  }
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

