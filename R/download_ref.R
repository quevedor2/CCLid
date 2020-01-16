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
                            refFileName = NULL, verbose = TRUE) {
  ccl.table <- availableRefCCL(saveDir = saveDir)
  whichx <- match(name, pSetTable[, 1])
  if (is.na(whichx)) {
    stop("Unknown Dataset. Please use the availableRefCCL() function for the table of available CCLid reference datasets.")
  }
  if (!file.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  if (is.null(refFileName)) {
    refFileName <- paste0(ccl.table[whichx, "Datasets"], ".rds")
  }
  if (!file.exists(file.path(saveDir, refFileName))) {
    downloader::download(url = as.character(ccl.table[whichx, "URL"]), 
                         destfile = file.path(saveDir, refFileName), 
                         quiet = !verbose)
  }
  ref.mat <- readRDS(file.path(saveDir, refFileName))
  return(get(ref.mat))
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
availableRefCCL <- function (saveDir = file.path(".", "CCLid"), 
                             myfn = "downloadTable.csv", verbose = TRUE) {
  if (!file.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  saveDir <- file.path("~", "git", "CCLid", "data-raw")
  dl.table <- read.csv(file.path(saveDir, myfn), check.names = FALSE, stringsAsFactors = FALSE)
  return(dl.table)
}