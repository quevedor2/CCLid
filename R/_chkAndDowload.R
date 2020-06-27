#' Check and download
#' @description Checks for an existing .rda file of the input
#' data type. If it doesn't exist, it downloads from the 
#' given CCLid::ccl_table data table URL and loads the .rda
#' file into memory
#' @param name Name of datatype to download and save to
#' @param whichx Row index of datatype to dowload
#' @param saveDir Directory that contains or will contain the .rda file
#' @return
#' Loads the .rda into .GlobalEnv
.chkAndDownload <- function(name, whichx, saveDir){
  rda_file <- file.path(saveDir, paste0(tolower(name), ".rda"))
  if(!file.exists(rda_file)){
    message(paste0("Downloading ", basename(rda_file), "..."))
    downloader::download(url = as.character(CCLid::ccl_table[whichx, "URL"]), 
                         destfile = rda_file, 
                         quiet = !verbose)
  } 
  load(rda_file, envir = .GlobalEnv)
}