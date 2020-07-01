#' getCVCL
#' @description get the Cellosaurus CVCL_ id for a given
#' cell line name
#'
#' @param cellid character: cell line name (e.g. MCF-7)
#' @param prioritize.datasets default=TRUE
#' @param melt.cells Melted cellosaurus dataframe, accessible from CCLid::ccl_table
#' 
#' @importFrom utils data
#'
#' @return CVCL_ style character name
#' @export
getCVCL <- function(cellid, prioritize.datasets=TRUE, melt.cells){
  # cellid <- "ES-2"
  #data(melt.cells)
  cl.match <- melt.cells[grep(paste0("^", cellid, "$"), melt.cells$ID),]

  ## If all CVCL's are equal, just use it
  if(length(unique(cl.match$CVCL))==1) cl.match <- cl.match[1,]


  if(nrow(cl.match) > 1 & prioritize.datasets) {
    ## Check if there is a single column that matches priority list
    priority.check <- colnames(cl.match)[9:ncol(cl.match)]
    row.s <- rowSums(cl.match[,priority.check])
    if(any(row.s > 0) & sum(as.logical(row.s)) == 1){
      cl.match <- cl.match[which(row.s > 0),]
    }
  }

  if(nrow(cl.match) > 1) {
    ## If no conclusion reached, throw everything back
    warning("Multiple cells found for the given ID (likely contamination), please manully select 1 CVCL for further analysis")
    return(cl.match)
  } else {
    return(cl.match$CVCL)
  }
}
