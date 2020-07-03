#' sigDiffBaf
#' @returns GRanges object of significant BAF drift (z > 3)
#' for a cna.obj returned by bafDrift$cna.obj
#' 
#' @param each.sample Each sample
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' 
#' @export
sigDiffBaf <- function(each.sample){
  es <- each.sample$output
  sig.idx <- which(es$t >= 3)
  sig.es=NULL
  if(length(sig.idx) > 0){
    es <- es[sig.idx,] 
    sig.es <- lapply(split(es, es$ID), makeGRangesFromDataFrame, keep.extra.columns=TRUE)
  }
  return(sig.es)
}