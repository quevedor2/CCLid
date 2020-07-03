#' splitConcordanceVals
#' @description Melts the data structure into Sample 1, Sample 2, concordance
#' and separates into Match/nonmatch if a reference metadata is given
#' 
#' @param dm dataframe melt
#' @param meta.df metadata for the melted dataframe
#'
#' @export
splitConcordanceVals <- function(dm, meta.df){
  dr.nm <- dm
  
  ## Changes matrix position [2,1] to a vector [2]
  .mPosToVec <- function(x, nrw){
    as.integer(sapply(nrw * (x-1), function(i) i+x))
  }
  
  n.rows <- nrow(dm)
  all.idx <- findCclPairs(meta.df, dr.nm)
  all.idx.v <- unlist(sapply(all.idx, .mPosToVec, nrw=n.rows))
  
  m.vals <- lapply(all.idx, function(idx){
    .meltDf(dr.nm[idx,idx,drop=FALSE])
  })
  dr.nm[all.idx.v] <- NA
  
  ## Old, slightly slower method:
  # m.vals <- list()
  # for(i in meta.df$ID){
  #   idx <- grep(paste0("_", i, "$"), colnames(dr.nm))
  #   m.val <- .meltDf(dr.nm[idx,idx,drop=FALSE])
  #   dr.nm[idx,idx] <- NA
  #   m.vals[[i]] <- m.val
  # }
  m.vals <- do.call(rbind, m.vals)
  nm.vals <- .meltDf(dr.nm)
  
  return(list("M"=m.vals,
              "NM"=nm.vals))
}