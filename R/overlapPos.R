#' .overlapProbeset
#' @description Overlaps probesets and orders based on chr position
#' 
#' @param ref.ids Reference probeset IDs to match order to
#' @param comp.ids IDs of the given input, comparator
#' 
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
#' @param comp Input comparison data with a $Probe_Set_ID column
#' @param ref Reference probeset order where rownames are Probe_Set_IDs
#' @param mapping If "probeset", will using $Probe_Set_ID from comp dataframe.
#' Else, it will use the rownames of ref as the probe set IDs
#'
#' @export
overlapPos <- function(comp, ref, mapping='probeset'){
  switch(mapping,
         "probeset"={
           if(any(grepl("Probe_Set_ID", colnames(comp)))){
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
#' @param meta.df Metadata containing cell IDs and their dataset mapping
#' @param dr.nm  matrix where column names are cell lines/datasets
#' @param ds Dataset (if any) to return (e.g. c('CCLE', 'GDSC'))
#' @importFrom utils setNames
#' 
#' @return A list of indices indicating pair of cell lines for the data matrix
#' @export
findCclPairs <- function(meta.df, dr.nm, ds=NULL){
  all.idx <- sapply(meta.df$ID, function(i) {
    id <- grep(paste0("_", i, "$"), x=colnames(dr.nm))
    id <- setNames(id, colnames(dr.nm)[id])
    if(!is.null(ds)){
      #ds <- c('CCLE', 'GDSC')
      id[grep(paste(ds, collapse="|"), names(id))]
    } else {
      id
    }
  })
  return(all.idx)
}
