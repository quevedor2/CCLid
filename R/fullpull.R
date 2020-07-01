#' fullpull
#' @description Gets the synonymous cell lines and cell line that they derived from,
#' as well as any indication of MSI or probelmatic cell line (e.g. contamination)
#'
#' @param cvcl Cellosaurus CVCL id (e.g. CVCL_1384)
#' @param melt.cells Melted cellosaurus dataframe, accessible from CCLid::ccl_table
#'
#' @return Dataframe of all synonymous (SS), similar origins (OI), MSI positive, or problematic
#' lines associated with the CVCL
#' @export
fullpull <- function(cvcl, melt.cells){
  if(!is.na(cvcl)){
    if(!substr(cvcl, 1, 5) == "CVCL_"){
      cvcl <- getCVCL(cvcl, melt.cells=melt.cells)
      if(class(cvcl) == 'data.frame' | length(cvcl) == 0) stop("Could not find the CVCL id for given input")
    }
  }
  
  tryCatch({
    rbind(.getSynonymous(cvcl, melt.cells=melt.cells),
          .getDerivedFrom(cvcl, melt.cells=melt.cells))
  }, error=function(e){NULL})
}
