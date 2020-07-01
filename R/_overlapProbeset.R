#' .overlapProbeset
#' @description Overlaps probesets and orders based on chr position
#' 
#' @param ref.ids Reference probeset IDs to match order to
#' @param comp.ids IDs of the given input, comparator
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#' 
.overlapProbeset <- function(ref.ids, comp.ids, snp6.dat){
  probe.meta <- snp6.dat$SNP$Probe_Set_ID
  
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