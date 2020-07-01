#' pIDs
#' @description Creates a standardized set of group names for 
#' samples using a reference
#' 
#' @param mat Matrix of BAF/Geno
#' @param meta.df Metadata for the samples being used with CVCL ids
#' @param datasets Datasets to include in analysis (Default 'GNE', 'CCLE', 'GDSC')
#'
#' @export
assignGrpIDs <- function(mat, meta.df, datasets=c("GNE", "GDSC", "CCLE")){
  # cvcl.idx <- grep("CVCL", colnames(meta.df))
  ds_idx <- which(colnames(meta.df) %in% c(datasets))
  meta_ds <- meta.df[,ds_idx]
  
  new.ids <- sapply(colnames(mat), function(i){
    i <- gsub("\\.[xy]$", "", i)
    suffix <- c("", ".cel", ".Cel", ".CEL")
    m_stat <- FALSE
    suffix_idx <- 1
    
    while(!m_stat){
      id_match <- which(paste0(i, suffix[suffix_idx])==meta_ds, arr.ind=TRUE)
      if(nrow(id_match) > 0) m_stat <- TRUE else  suffix_idx <- suffix_idx+1
      if(suffix_idx > 4) m_stat <- TRUE
      #print(paste(suffix[suffix_idx], id_match))
    }
    
    if(nrow(id_match) > 0){
      ridx <- id_match[1,1]
      cidx <- id_match[1,2]
      paste0(colnames(meta_ds)[cidx], "_", meta.df[ridx,]$ID)
    } else {
      i
    }
  })
  return(new.ids)
}
