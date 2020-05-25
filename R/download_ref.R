#' downloadRefCCL
#' @description Downloads the precomputed RDS data structures for reference
#' 
#' @param name 
#' @param saveDir 
#' @param verbose 
#' @param refFileName 
#'
#' @return
#' @export
#'
#' @examples
#' downloadRefCCL(name="BAF")
downloadRefCCL <- function (name, saveDir = file.path(".", "CCLid"), 
                            refFileName = NULL, verbose = FALSE, bin.size=NULL) {
  ccl.table <- CCLid::availableRefCCL(saveDir = saveDir)
  whichx <- match(name, ccl.table[, 1])
  if (is.na(whichx)) {
    stop("Unknown Dataset. Please use the availableRefCCL() function for the table of available CCLid reference datasets.")
  }
  if (!file.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  if (is.null(refFileName)) {
    refFileName <- paste0(ccl.table[whichx, "Ref.type"], ".rds")
  }

  if(file.exists(file.path(saveDir, paste0("ref_", as.integer(bin.size), ".desc")))){
    ## Looks for pre-existing bigmemory data structure to circumvent loading into memory
    require(bigmemory)
    require(biganalytics)
    
    if(verbose) print("Reading in existing bigmemory object...")
    shared.desc <- dget(file.path(saveDir, paste0("ref_", as.integer(bin.size), ".desc")))
    shared.bigobject <- attach.big.matrix(shared.desc)
    ids <- readRDS(file.path(saveDir, "ref_mat_ID.rds"))
    ref.mat <- shared.bigobject
    options(bigmemory.allow.dimnames=TRUE)
    rownames(ref.mat) <- ids
  } else {
    if (!file.exists(file.path(saveDir, refFileName))) {
      downloader::download(url = as.character(ccl.table[whichx, "URL"]), 
                           destfile = file.path(saveDir, paste0(refFileName, ".gz")), 
                           quiet = !verbose)
      if(verbose) print(paste0("Unzipping: ", file.path(saveDir, paste0(refFileName, ".gz"))))
      system(command = paste0('gunzip ', file.path(saveDir, paste0(refFileName, ".gz"))))
    }
    
    ref.mat <- readRDS(file.path(saveDir, refFileName))
  }
  
  return(ref.mat)
}

#' availableRefCCL
#' @description table describing the different datasets for download 
#' 
#' @param saveDir 
#' @param myfn 
#' @param verbose 
#'
#' @return
#' @export
availableRefCCL <- function (saveDir = file.path(".", "CCLid"), tableDir=NULL,
                             myfn = "downloadTable.csv", verbose = TRUE) {
  if(is.null(tableDir)){
    tableDir <- system.file(file.path("extdata"), package="CCLid")   #file.path("~", "git", "CCLid", "data-raw")
  }
  if (!file.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  dl.table <- read.csv(file.path(tableDir, myfn), check.names = FALSE, stringsAsFactors = FALSE)
  return(dl.table)
}

#' Title
#'
#' @param ref.mat 
#' @param analysis 
#'
#' @return
#' @export
#'
#' @examples
formatRefMat <- function(name, ref.mat, analysis, 
                         varFileName=NULL, saveDir = file.path(".", "CCLid"), 
                         bin.size=1e6, just.var=FALSE, fill.na=FALSE, verbose=FALSE){
  # saveDir=PDIR
  # name='BAF'
  # analysis='baf'
  # just.var=FALSE
  # fill.na=FALSE
  # varFileName=NULL
  ## Set filename and file path
  ccl.table <- CCLid::availableRefCCL(saveDir = saveDir)
  whichx <- match(name, ccl.table[, 1])
  
  if (is.null(varFileName)) {
    varFileName <- paste0(as.integer(bin.size), ".", ccl.table[whichx, "Ref.type"], ".rds")
  }
  
  ## Process if existing RDS doesn't exist
  if(!just.var){
    if(any(grepl("^ID$", colnames(ref.mat)))) {
      rownames(ref.mat) <- ref.mat$ID
      ref.mat <- ref.mat[,-1]
    }
    keep.idx <- switch(analysis,
                       lrr=grep("CN", gsub("_.*", "", rownames(ref.mat))),
                       baf=grep("SNP", gsub("_.*", "", rownames(ref.mat))),
                       stop("'analysis' parameter must be submitted: 'lrr, baf'"))
    ref.mat <- ref.mat[keep.idx,]
    if(fill.na) ref.mat[is.na(ref.mat)] <- median(as.matrix(ref.mat), na.rm=T)
  }
  
  ## Calculate variant features if file doesn't already exist
  if (!file.exists(file.path(saveDir, varFileName))) {
    if(verbose) print("Generate feature variance data")
    var.feats <- .getVariantFeatures(ref.mat, bin.size)
    saveRDS(var.feats, file.path(saveDir, varFileName))
  } else {
    if(verbose) print("Reading existing variance data")
    var.feats <- readRDS(file.path(saveDir, varFileName))
  }
  
  return(list("mat"=ref.mat,
              "var"=var.feats))
}

#' .getVariantFeatures
#' @description Calculates the features/probesets with the most variance from a 
#' given matrix (rownames).  Returns equally spaced Variant Probesets
#' 
#' @param ref.mat 
#' @param bin.size 
#'
#' @return
.getVariantFeatures <- function(ref.mat, bin.size=1e6){
  require(BSgenome.Hsapiens.UCSC.hg19)
  
  ## Gets variance of matrix
  ref.vars <- apply(ref.mat, 1, var, na.rm=TRUE)
  ref.vars <- setNames(round(ref.vars,3), rownames(ref.mat))
  
  ## Get tally of non-NA sites:
  na.per.snp <- apply(ref.mat, 1, function(i) sum(is.na(i)))
  na.per.snp <- setNames(na.per.snp, rownames(ref.mat))
  
  ## Identifies genomic position of SNPs in ref.mat
  all.pb <-  as.data.frame(CCLid::snp6.dat$All)[,c(1:3)] #Loads in probeset meta data
  rownames(all.pb) <- CCLid::snp6.dat$All$Probe_Set_ID
  var.pb <- all.pb[names(ref.vars),] # Selects rows based on probeset IDs
  var.pb$var <- round(ref.vars,3)  # Variance column
  var.pb$Probe_Set_ID <- names(ref.vars)  # Probe_Set_ID column
  var.pb$num.snps <- na.per.snp # Number of samples that have BAF value
  na.idx <- which(is.na(var.pb$seqnames)) #Removes NA if any
  if(length(na.idx) > 0) var.pb <- var.pb[-na.idx,]
  var.gr <- sort(makeGRangesFromDataFrame(var.pb, keep.extra.columns = TRUE))
  
  ## Cuts the chromosomes into equal sized bins
  chr.sizes <- seqlengths(Hsapiens)[paste0("chr", c(1:22, "X", "Y"))]
  bins   <- tileGenome(chr.sizes, tilewidth=bin.size, cut.last.tile.in.chrom=T)
  seqlevelsStyle(bins) <- seqlevelsStyle(var.gr) <- 'UCSC'
  
  ## Finds the overlap of SNP probesets with the genomic bins
  ov.idx <- findOverlaps(bins, var.gr)
  ov.split <- split(ov.idx, queryHits(ov.idx))
  
  ## For each bin, selects the Probeset with the max variance and returns that
  var.l <- lapply(ov.split, function(i){
    binned.var <- var.gr[subjectHits(i),]
    n.var <- as.data.frame(binned.var)[,c("var", "num.snps")]
    n.var
    # max.idx <- which.max(search.space$var)
    # search.space[max.idx,]
  })
  
  return(var.l)
}

#' pIDs
#' @description Creates a standardized set of group names for 
#' samples using a reference
#' 
#' @param mat 
#' @param meta.df 
#'
#' @return
#' @export
assignGrpIDs <- function(mat, meta.df){
  cvcl.idx <- grep("CVCL", colnames(meta.df))
  new.ids <- sapply(colnames(mat), function(i){
    i <- gsub("\\.[xy]$", "", i)
    cidx <- grep(i, meta.df[,-cvcl.idx], ignore.case = T)
    cidx <- cidx[length(cidx)]
    ridx <- grep(paste0("^", i, "(.cel)?$"), meta.df[,cidx], ignore.case = T)
    
    if(length(ridx) > 0){
      paste0(colnames(meta.df)[cidx], "_", meta.df[ridx,]$ID)
    } else {
      i
    }
  })
  return(new.ids)
}