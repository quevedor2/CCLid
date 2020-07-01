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
    print(head(meta.df))
    print("Assigning group IDs")
    new.ids <- assignGrpIDs(ref.mat, meta.df)
    new.ids[duplicated(new.ids)] <- gsub("_", "2_",new.ids[duplicated(new.ids)])
    save(new.ids, file=file.path(PDIR, "col_ids.rda"))
  } else {
    print("Loading in existing group IDs")
    load(file=file.path(PDIR, "col_ids.rda"))
  }
  colnames(ref.mat) <- new.ids
  return(list("ref"=ref.mat,
              "var"=var.dat))
}

#' compareVcf
#' @description Checks a VCF file(s) against reference dataset
#' to look for similarity to any known cell lines
#'
#' @param vcfFile Path to VCF file to check against reference datasets
#' @param var.dat Variance data (list)
#' @param ref.mat Reference matrix (matrix)
#' @param max.snps Max number of SNPs to reduce 
#' @param ids IDs
#' @param sampletype Strictly for labelling purposes
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#' @param ... Extra param
#' @param all.ids IDs to subset to
#'
#' @return Matrix: Containing reference matrix subsetted to common SNPs as 
#' the input VCF, as well as a left-joined VCF data
#' @export
#'
compareVcf <- function(vcfFile, var.dat, ref.mat, 
                       max.snps=1e6, ids=NULL, sampletype='RNA', all.ids=NULL, 
                       snp6.dat, ...){
  vcf.map <- CCLid::mapVcf2Affy(vcfFile, snp6.dat=snp6.dat)
  vcf.map <- .filt(vcf.map, ...) ## Memory: up to 1.8Gb 
  
  ## Combine matrices and reduce features
  ## Find the overlap between the COMParator and the REFerence
  vcf.map.var <- CCLid::mapVariantFeat(vcf.map, var.dat)
  vcf.to.use <- vcf.map.var
  ov.idx <- CCLid::overlapPos(comp = vcf.to.use$BAF,
                              ref=ref.mat, mapping = 'probeset')
  if(nrow(ov.idx) > max.snps){
    ov.idx <- ov.idx[order(ov.idx$ref)[1:max.snps],]
  }
  rm(vcf.map, vcf.map.var); gc()  ## Cleanup
  
  ## BIG memory sink: shoots up to 7Gb
  if(is.null(ids)){
    refm <- ref.mat
  } else {
    colidx <- which(colnames(ref.mat) %in% all.ids)
    message("Isolating for: ", paste(colnames(ref.mat)[colidx], collapse=", "))
    refm <- ref.mat[,colidx,drop=FALSE]
  }
  x.mat <- cbind(vcf.to.use$BAF$BAF[ov.idx$comp], 
                 refm[ov.idx$ref,])
  colnames(x.mat)[1] <- paste0(sampletype, "_", gsub(".vcf.*", "", basename(vcfFile)))
  
  if(storage.mode(ref.mat[,1]) == 'integer'){
    x.mat[,-1] <- x.mat[,-1] / 100
  }
  return(x.mat)
}

