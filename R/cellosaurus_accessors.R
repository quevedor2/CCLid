## Private Functions
# accessor for OI/SS

#' .getDerivedFrom
#' @param cvcl CVCL id
#'
.getDerivedFrom <- function(cvcl){
  cvcl.oi <- unique(melt.cells[grep(paste0("^", cvcl, "$"), melt.cells$CVCL),]$OI)
  oi.list <- lapply(cvcl.oi, function(cv) melt.cells[grep(paste0("^", cv, "$"), melt.cells$CVCL),])
  cl.match <- do.call(rbind, oi.list)
  if(nrow(cl.match) > 0){
    cl.match$acr <- "OI"
  }
  return(cl.match)
}

#' .getSynonymous
#' @param cvcl CVCL id
#'
.getSynonymous <- function(cvcl){
  cl.match <- melt.cells[grep(paste0("^", cvcl, "$"), melt.cells$CVCL),]
  if(nrow(cl.match) > 0){
    cl.match$acr <- "SS"
  }
  return(cl.match)
}
