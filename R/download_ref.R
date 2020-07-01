#' .getVariantFeatures
#' @description Calculates the features/probesets with the most variance from a 
#' given matrix (rownames).  Returns equally spaced Variant Probesets
#' 
#' @param ref.mat Reference matrix of Samples by Probesets
#' @param bin.size Bin size, Default=1e6
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
.getVariantFeatures <- function(ref.mat, bin.size=1e6){
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

#' pIDs
#' @description Creates a standardized set of group names for 
#' samples using a reference
#' 
#' @param mat Matrix of BAF/Geno
#' @param meta.df Metadata for the samples being used with CVCL ids
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
