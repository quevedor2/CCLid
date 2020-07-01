#' Check and download
#' @description Checks for an existing .rda file of the input
#' data type. If it doesn't exist, it downloads from the
#' given CCLid::ccl_table data table URL and loads the .rda
#' file into memory
#' @param name Name of datatype to download and save to
#' @param whichx Row index of datatype to dowload
#' @param saveDir Directory that contains or will contain the .rda file
#' @param verbose Verbose
#' @return
#' Loads the .rda into .GlobalEnv
#' @export
#' @keywords internal
.chkAndDownload <- function(name, whichx, saveDir,  verbose=FALSE){
  rda_file <- file.path(saveDir, paste0(tolower(name), ".rda"))
  if(!file.exists(rda_file)){
    message(paste0("Downloading ", basename(rda_file), "..."))
    downloader::download(url = as.character(CCLid::ccl_table[whichx, "URL"]),
                         destfile = rda_file,
                         quiet = !verbose)
  }
  load(rda_file, envir = .GlobalEnv)
  #load(rda_file, envir=env)
  #return(env)
}

#' .genMatchSnpDiff
#' @description generates a matrix of euclidean distance between probesets of 
#' sample pairs if the sample id's MATCH
#'
#' @param meta.df Cell line metadata, accessible from CCLid::ccl_table
#' @param dr.nm nonmatch SNP
#'
#' @importFrom utils combn
#' @importFrom matrixStats rowDiffs
#' @export
#' @keywords internal
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

#' .getVariantFeatures
#' @description Calculates the features/probesets with the most variance from a 
#' given matrix (rownames).  Returns equally spaced Variant Probesets
#' 
#' @param ref.mat Reference matrix of Samples by Probesets
#' @param bin.size Bin size, Default=1e6
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#' 
#' @importFrom stats var
#' @importFrom stats setNames
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @importFrom GenomicRanges tileGenome
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @export
#' @keywords internal
.getVariantFeatures <- function(ref.mat, bin.size=1e6, snp6.dat){
  ## Gets variance of matrix
  ref.vars <- apply(ref.mat, 1, var, na.rm=TRUE)
  ref.vars <- setNames(round(ref.vars,3), rownames(ref.mat))
  
  ## Get tally of non-NA sites:
  na.per.snp <- apply(ref.mat, 1, function(i) sum(is.na(i)))
  na.per.snp <- setNames(na.per.snp, rownames(ref.mat))
  
  ## Identifies genomic position of SNPs in ref.mat
  all.pb <-  as.data.frame(snp6.dat$All)[,c(1:3)] #Loads in probeset meta data
  rownames(all.pb) <- snp6.dat$All$Probe_Set_ID
  var.pb <- all.pb[names(ref.vars),] # Selects rows based on probeset IDs
  var.pb$var <- round(ref.vars,3)  # Variance column
  var.pb$Probe_Set_ID <- names(ref.vars)  # Probe_Set_ID column
  var.pb$num.snps <- na.per.snp # Number of samples that have BAF value
  na.idx <- which(is.na(var.pb$seqnames)) #Removes NA if any
  if(length(na.idx) > 0) var.pb <- var.pb[-na.idx,]
  var.gr <- sort(makeGRangesFromDataFrame(var.pb, keep.extra.columns = TRUE))
  
  ## Cuts the chromosomes into equal sized bins
  chr.sizes <- seqlengths(Hsapiens)[paste0("chr", c(1:22, "X", "Y"))]
  bins   <- tileGenome(chr.sizes, tilewidth=bin.size, cut.last.tile.in.chrom=T)
  seqlevelsStyle(bins) <- seqlevelsStyle(var.gr) <- 'UCSC'
  
  ## Finds the overlap of SNP probesets with the genomic bins
  ov.idx <- findOverlaps(bins, var.gr)
  ov.split <- split(ov.idx, queryHits(ov.idx))
  
  ## For each bin, selects the Probeset with the max variance and returns that
  var.l <- lapply(ov.split, function(i){
    binned.var <- var.gr[subjectHits(i),]
    n.var <- as.data.frame(binned.var)[,c("var", "num.snps")]
    n.var
    # max.idx <- which.max(search.space$var)
    # search.space[max.idx,]
  })
  
  return(var.l)
}

#' .overlapProbeset
#' @description Overlaps probesets and orders based on chr position
#' 
#' @param ref.ids Reference probeset IDs to match order to
#' @param comp.ids IDs of the given input, comparator
#' @param snp6.dat SNP6 probeset genomic position, accessible from CCLid::ccl_table
#' 
#' @export
#' @keywords internal
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

#' .getDerivedFrom
#' @param cvcl CVCL id
#' @param melt.cells Melted cellosaurus dataframe, accessible from CCLid::ccl_table
#'
#' @export
#' @keywords internal
.getDerivedFrom <- function(cvcl, melt.cells){
  cvcl.oi <- unique(melt.cells[grep(paste0("^", cvcl, "$"), melt.cells$CVCL),]$OI)
  oi.list <- lapply(cvcl.oi, function(cv) melt.cells[grep(paste0("^", cv, "$"), melt.cells$CVCL),])
  cl.match <- do.call(rbind, oi.list)
  if(nrow(cl.match) > 0){
    cl.match$acr <- "OI"
  }
  return(cl.match)
}

#' .getSynonymous
#' @param cvcl CVCL id
#' @param melt.cells Melted cellosaurus dataframe, accessible from CCLid::ccl_table
#'
#' @export
#' @keywords internal
.getSynonymous <- function(cvcl, melt.cells){
  cl.match <- melt.cells[grep(paste0("^", cvcl, "$"), melt.cells$CVCL),]
  if(nrow(cl.match) > 0){
    cl.match$acr <- "SS"
  }
  return(cl.match)
}
