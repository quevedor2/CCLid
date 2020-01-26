#' .overlapProbeset
#' @description Overlaps probesets and orders based on chr position
#' 
#' @param ref.ids 
#' @param comp.ids 
#'
#' @return
.overlapProbeset <- function(ref.ids, comp.ids){
  probe.meta <- CCLid::snp6.dat$SNP$Probe_Set_ID
  
  # Order according the meta data for probesets
  idx.df <- data.frame("comp"=match(probe.meta, comp.ids),
                       "ref"=match(probe.meta, ref.ids))
  # Identify instances where probesets are found in both Ref and Comp
  non.na <- which(rowSums(is.na(idx.df)) == 0)
  
  if(length(non.na) > 0){
    idx.df[non.na,]
  } else {
    stop("Could not find any overlapping SNPs between datasets")
  }
}



#' overlapPos
#' @description parent function that handles how overlaps are done between COMParing 
#' group and REFerence datasets
#' 
#' @param comp 
#' @param ref 
#' @param mapping 
#'
#' @return
#' @export
overlapPos <- function(comp, ref, mapping='probeset'){
  switch(mapping,
         "probeset"={
           if(grepl("Probe_Set_ID", colnames(comp))){
             .overlapProbeset(ref.ids=rownames(ref), 
                                     comp.ids=comp$Probe_Set_ID)
           } else {
             .overlapProbeset(ref.ids=rownames(ref), 
                              comp.ids=rownames(comp))
           }
           })
}

#' findCclPairs
#' @description Finds indices for cell lines that are found
#' in multiple datasets
#' 
#' @param meta.df 
#' @param dr.nm  matrix where column names are cell lines/datasets
#'
#' @return
#' @export
findCclPairs <- function(meta.df, dr.nm){
  all.idx <- sapply(meta.df$ID, function(i) grep(paste0("_", i, "$"), x=colnames(dr.nm)))
  return(all.idx)
}