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



overlapPos <- function(comp, ref, mapping='probeset'){
  switch(mapping,
         "probeset"=.overlapProbeset(ref.ids=rownames(ref), 
                                     comp.ids=comp$Probe_Set_ID))
}
