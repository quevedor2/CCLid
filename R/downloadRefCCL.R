#' downloadRefCCL
#' @description Downloads the precomputed RDS data structures for reference
#' 
#' @param name Currently only supports "BAF"
#' @param saveDir Directory containign bigmemory .desc/.bin folder, or where to save new download
#' @param verbose Verbose (Default=FALSE)
#' @param ... Extra snp6.dat
#' @importFrom bigmemory attach.big.matrix
#' 
#' @export
#'
downloadRefCCL <- function (name, saveDir = file.path(".", "CCLid"), 
                            verbose = FALSE, ...) {
  if (!file.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  
  if(!any(toupper(name) %in% CCLid::ccl_table$`Data.type`)){
    stop("'name' must be one of the following data types: ", 
         paste(CCLid::ccl_table$`Data.type`, collapse=", "), ".\n", 
         "Use ?ccl_table for more information on the different data types.")
  }
  
  whichx <- match(toupper(name), CCLid::ccl_table$`Data.type`)
  refFileName <- paste0(toupper(name), ".rds")
  
  ## Load in the probeset BAF/Genotype matrix
  if(name %in% c("BAF", "GENO")){
    identifier <- paste0("ref_", tolower(name))
    row_ids <- paste0("ref_", tolower(name), "-ids.rds")
    rds_chk <- file.exists(file.path(saveDir, paste0(identifier, ".rds")))
    desc_chk <- file.exists(file.path(saveDir, paste0(identifier, ".desc")))
    bin_chk <- file.exists(file.path(saveDir, paste0(identifier, ".bin")))
    row_chk <- file.exists(file.path(saveDir, row_ids))
    
    ## Check if the files exist, else download them
    if (!(desc_chk & bin_chk & row_chk)) {
      message("Downloading .bin, .desc, and row_ids file...")
      # URL: .bin file
      if(!bin_chk){
        downloader::download(url = as.character(CCLid::ccl_table[whichx, "URL"]), 
                             destfile = file.path(saveDir, paste0(identifier, ".bin")), 
                             quiet = !verbose)
      }
      
      # URL2: .desc file
      if(!desc_chk){
        downloader::download(url = as.character(CCLid::ccl_table[whichx, "URL2"]), 
                             destfile = file.path(saveDir, paste0(identifier, ".desc")), 
                             quiet = !verbose)
      }
      
      # URL3: row IDs file
      if(!row_chk){
        downloader::download(url = as.character(CCLid::ccl_table[whichx, "URL3"]), 
                             destfile = file.path(saveDir, row_ids), 
                             quiet = !verbose)
      }
    }
    
    ## Read in bigmemory stuff     
    if(rds_chk){
      if(verbose) print("Reading in RDS file for sample x probeset object...")
      ref.mat <- readRDS(file.path(saveDir, paste0(identifier, ".rds")))
    } else {
      if(verbose) print("Reading in bigmemory sample x probeset object...")
      shared.desc <- dget(file.path(saveDir, paste0(identifier, ".desc")))
      shared.desc@description$dirname <- gsub("\\/?$", "/", saveDir)
      shared.desc@description$filename <- paste0(identifier, ".bin")
      shared.bigobject <- attach.big.matrix(shared.desc)
      ref.mat <- shared.bigobject
      options(bigmemory.allow.dimnames=TRUE)
    }
    
    
    ids <- readRDS(file.path(saveDir, row_ids))
    rownames(ref.mat) <- ids
    
    return(ref.mat)
  } else {
    .chkAndDownload(name, whichx, saveDir)
  }
}
