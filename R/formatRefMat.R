#' formatRefMat
#' @description Pre-processes the input data into a matrix of BAF values and 
#' the least/most variant SNPs within that cohort
#' 
#' @param ref.mat A reference matrix of all pharmacogenomic reference samples and their BAFs,
#' typically obtained from the downloadRefCCL() function
#' @param analysis Only 'baf' is implemented at this point
#' @param name Name to be applied to the returned object
#' @param varFileName RDS containing the variant SNP information
#' @param saveDir Directory to save/load the variant RDS file to
#' @param bin.size Default is set to 5e5, a variant file must be created for this bin size
#' @param just.var Just runs the variant SNP part of the script, skips subsetting
#' @param fill.na Fills NA with median (Default=FALSE)
#' @param verbose Verbose
#' @importFrom stats median
#' 
#' @return Returns a list object:
#' 'ref' = matrix of SNPs by samples for least variant SNPs
#' 'var' = list of each 'bin size' and the SNPs and their BAF that populate it
#' @export
#'
formatRefMat <- function(name, ref.mat, analysis='baf', 
                         varFileName=NULL, saveDir = file.path(".", "CCLid"), 
                         bin.size=1e6, just.var=FALSE, fill.na=FALSE, verbose=FALSE){
  # saveDir=PDIR
  # name='BAF'
  # analysis='baf'
  # just.var=FALSE
  # fill.na=FALSE
  # varFileName=NULL
  ## Set filename and file path
  whichx <- match(toupper(name), CCLid::ccl_table[, 1])
  
  if (is.null(varFileName)) {
    varFileName <- paste0(as.integer(bin.size), ".", 
                          CCLid::ccl_table[whichx, "Ref.type"], ".rds")
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
    which_file <- match(varFileName, CCLid::ccl_table[, 1])
    if(!is.na(which_file)){
      print("Downloading and loading existing variance file...")
      downloader::download(url = as.character(CCLid::ccl_table[which_file, "URL"]), 
                           destfile = file.path(saveDir, varFileName), 
                           quiet = !verbose)
      var.feats <- readRDS(file.path(saveDir, varFileName))
    } else {
      if(verbose) print("Generate feature variance data")
      var.feats <- .getVariantFeatures(ref.mat, bin.size)
      saveRDS(var.feats, file.path(saveDir, varFileName))
    }
  } else {
    if(verbose) print("Reading existing variance data")
    var.feats <- readRDS(file.path(saveDir, varFileName))
  }
  
  return(list("mat"=ref.mat,
              "var"=var.feats))
}