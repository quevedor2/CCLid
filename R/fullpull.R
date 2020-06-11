#' fullpull
#' @description Gets the synonymous cell lines and cell line that they derived from,
#' as well as any indication of MSI or probelmatic cell line (e.g. contamination)
#'
#' @param cvcl Cellosaurus CVCL id (e.g. CVCL_1384)
#' @param melt.cells data.frame made from cellosaurus xml
#'
#' @return Dataframe of all synonymous (SS), similar origins (OI), MSI positive, or problematic
#' lines associated with the CVCL
#' @export
#'
#' @examples data(melt.cells)
#' cvcl <- getCVCL('Hela', melt.cells)
#' fullpull(cvcl, melt.cells)
#'
#' fullpull("HeLa", melt.cells)
fullpull <- function(cvcl, melt.cells){
  if(!is.na(cvcl)){
    if(!substr(cvcl, 1, 5) == "CVCL_"){
      cvcl <- getCVCL(cvcl, melt.cells)
      if(class(cvcl) == 'data.frame' | length(cvcl) == 0) stop("Could not find the CVCL id for given input")
    }
  }
  
  tryCatch({
    rbind(.getSynonymous(cvcl, melt.cells),
          .getDerivedFrom(cvcl, melt.cells))
  }, error=function(e){NULL})
}
