#' loadRef
#' @description Wrapper to load in all the reference matrices
#' and annotate them with groups, as well as calculate variance
#'
#' @param analysis Only support BAF at the moment
#' @param ... Extra param
#' @param PDIR Directory for download data
#' @param rm.gne Remove gCSI from the analysis (default=FALSE)
#' @param bin.size Bin size (default=1e6)
#' @param meta.df Cell line metadata, accessible from CCLid::ccl_table
#' @param verbose Verbose
#' 
#' @return List containing "ref"=reference BAF matrix
#' and "var"=list containing variance for bin-sizes
#' @export
#'
loadRef <- function(PDIR=NULL, analysis='baf', rm.gne=FALSE, 
                    bin.size=1e6, verbose=FALSE, meta.df=NULL, ...){
  # if(verbose) print("Checking for existing reference data...")
  # ref.data.exists <- any(grepl(paste0(as.integer(bin.size), ".", toupper(analysis)), list.files(PDIR)))
  # 
  # if(ref.data.exists){
  #   if(verbose) print(paste0("Reading existing data: ", paste0(as.integer(bin.size), ".", toupper(analysis), ".rds"), "..."))
  #   ref.dat <- file.path(PDIR, paste0(as.integer(bin.size), ".", toupper(analysis), ".rds"))
  #   format.dat <- readRDS(ref.dat)
  # } else {
  # 
  # }
  #if(verbose) print("Downloading and loading metadata...")
  #metadata <- c("meta.df",  "affy.omni", "cin70", "gne.meta", "melt.cells", "snp6.dat")
  #sapply(metadata, downloadRefCCL, saveDir=PDIR, verbose=verbose)
  #env <- new.env()
  #for(m in metadata){
  #  downloadRefCCL(name=m, saveDir=PDIR, env=env, verbose=verbose)
  #}
  
  if(verbose) print("Attaching/downloading reference matrix...")
  ref.mat <- downloadRefCCL(toupper(analysis), saveDir = PDIR, verbose=verbose)
  
  if(rm.gne){
    rm.idx <- grep("^Unk*", colnames(ref.mat))
    ref.mat <- ref.mat[,-rm.idx]
  }
  
  if(verbose) print("Ensuring proper format and caluclating variance per bin...")
  format.dat <- formatRefMat(name=toupper(analysis), ref.mat=ref.mat, saveDir = PDIR, 
                             analysis=tolower(analysis), bin.size=bin.size, ...) #bin.size=5e5
  
  ref.mat <- format.dat$mat
  var.dat <- format.dat$var
  rm(format.dat)
  
  ## Assign group IDs (e.g. 22Rv1.cel -> GDSC_22Rv1)
  if(verbose) print("Assigning group IDs...")
  if(!file.exists(file=file.path(PDIR, "col_ids.rda"))){
    if(verbose) print(head(meta.df))
    if(verbose) print("Assigning group IDs")
    new.ids <- assignGrpIDs(ref.mat, meta.df)
    new.ids[duplicated(new.ids)] <- gsub("_", "2_",new.ids[duplicated(new.ids)])
    save(new.ids, file=file.path(PDIR, "col_ids.rda"))
  } else {
    if(verbose) print("Loading in existing group IDs")
    load(file=file.path(PDIR, "col_ids.rda"))
  }
  colnames(ref.mat) <- new.ids
  return(list("ref"=ref.mat,
              "var"=var.dat))
}