#' .genMatchSnpDiff
#' @description generates a matrix of euclidean distance between probesets of 
#' sample pairs if the sample id's MATCH
#'
#' @param meta.df Cell line metadata, accessible from CCLid::ccl_table
#' @param dr.nm nonmatch SNP
#'
#' @importFrom utils combn
#' @importFrom matrixStats rowDiffs
.genMatchSnpDiff <- function(dr.nm, meta.df){
  #data(meta.df)
  m.d <- sapply(meta.df$ID, function(i){
    idx <- grep(paste0("_", i, "$"), colnames(dr.nm))
    if(length(idx) > 1){
      d <- as.data.frame(apply(combn(idx,2), 2, function(com){
        rowDiffs(as.matrix(dr.nm[,com]))
      }))
      rownames(d) <- rownames(dr.nm)
      colnames(d) <- apply(combn(idx,2),2, function(j) paste(colnames(dr.nm)[j],collapse=':'))
      return(d)
    }
  })
  null.idx <- sapply(m.d, is.null)
  if(any(null.idx)){
    m.d <- do.call(cbind, m.d[-which(null.idx)])
  } 
  
  return(m.d)
}