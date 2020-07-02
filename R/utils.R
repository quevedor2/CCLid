#' .filt
#' @description Removes probesets under a certain depth, especially useful
#' for RNAsequencing vcfs
#'
#' @param vcf list from mapVcf2Affy
#' @param min.depth Default=5
#'
#' @return A filtered list from mapVcf2Affy
#' @export
#' @keywords internal
.filt <- function(vcf, min.depth=5){
  rm.idx <- which(vcf[[1]]$tdepth <= min.depth)
  vcf[[1]] <- vcf[[1]][-rm.idx,,drop=FALSE]
  vcf[[2]] <- vcf[[2]][-rm.idx,,drop=FALSE]
  return(vcf)
}

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
  ## Gets variance of matrix | Tally of NA SNPs
  if(class(ref.mat) == 'big.matrix'){
    ref.vars <- biganalytics::apply(ref.mat, 1, var, na.rm=TRUE)
    na.per.snp <- biganalytics::apply(ref.mat, 1, function(i) sum(is.na(i)))
  } else {
    ref.vars <- apply(ref.mat, 1, var, na.rm=TRUE)
    na.per.snp <- apply(ref.mat, 1, function(i) sum(is.na(i)))
  }
  
  ref.vars <- setNames(round(ref.vars,3), rownames(ref.mat))
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

#' .vcf2AffyGT
#' @description Maps Genotypes called from VCFs (0/0, 0/1, 1/0, 1/1) to the Affymetrix
#' standard of 0, 1, 2 genotypes
#' 
#' @param vcf.gt VariantAnnotation object where the $GT column contains 1/1, 0/1, etc...
#'
#' @export
#' @keywords internal
.vcf2AffyGT <- function(vcf.gt){
  vcf.affy.map <- c('1/1'=2, 
                    '0/1'=1,
                    '1/0'=1,
                    '0/0'=0)
  if(class(vcf.gt) == 'VRanges'){
    affy.gt <- vcf.affy.map[vcf.gt$GT]
  } else if (all(vcf.gt %in% names(vcf.gt))) {
    affy.gt <- vcf.affy.map[vcf.gt]
  } else {
    stop("Could not determine type of input data, must be vector of 0/0 genotypes or VariantAnnotation object")
  }
  if(any(is.na(affy.gt))) affy.gt[is.na(affy.gt)] <- -1
  
  return(affy.gt)
}

#' .revComp
#' @description Reverse complements the Alleles
#' 
#' @param vcf.gt VariantAnnotation ob ject
#' @param ret.idx  Boolean, if True, does not reverse complement any alleles
#' @export
#' @keywords internal
.revComp <- function(vcf.gt, ret.idx=FALSE){
  complementAllele <- c(A="T", T="A", C="G", G="C")
  neg.strand <- which(vcf.gt$Strand == '-')
  
  if(ret.idx){
    neg.strand
  } else {
    vcf.gt[neg.strand,]$Allele_A <- complementAllele[vcf.gt[neg.strand,]$Allele_A]
    vcf.gt[neg.strand,]$Allele_B <- complementAllele[vcf.gt[neg.strand,]$Allele_B]
    vcf.gt
  }
}

#' .fixGT
#' @description Switches the genotype to adjust for the reverse complement switch (0/0 -> 1/1)
#' 
#' @param vcf.gt VCF genotype (0,1,2)
#' @param affy.gt Affy6 genotype
#' @param ret.idx Returns the index of matching instead of genotype fix
#' @importFrom VariantAnnotation alt
#' @export
#' @keywords internal
.fixGT <- function(vcf.gt, affy.gt, ret.idx=FALSE){
  complementGenotype <- c('0'='2', 
                          '1'='1', 
                          '2'='0')
  rev.idx <- which(vcf.gt$Allele_A == alt(vcf.gt))
  if(ret.idx){
    rev.idx
  } else {
    affy.gt[rev.idx] <- complementGenotype[as.character(affy.gt[rev.idx])]
    as.integer(affy.gt)
  }
}

#' .normBAF
#' @description rescales the BAF from a 0-1 scale to a 0-0.5
#' 
#' @param x BAF vector
#' @param lower if TRUE, 0-0.5 range, ELSE, 0.5-1 range
#' @export
#' @keywords internal
.normBAF <- function(x, lower=T){
  #x <- seq(0, 1, by=0.1)
  x[x>1] <- 1; x[x<0] <- 0
  if(lower){
    nBAF <- (0.5 - abs(x - 0.5))
  } else {
    nBAF <- (0.5 + abs(x - 0.5))
  }
  
  nBAF
}

#' .jsonToGr
#' @description converts a JSON object into a VariantAnnotation object
#'
#' @param from.file Boolean to say whether it comes from a file or passed in object
#' @param json Json VCF format - OUTDATED
#' @importFrom VariantAnnotation VRanges
#' @importFrom IRanges IRanges
#' @export
#' @keywords internal
.jsonToGr <- function(json, from.file=TRUE){
  if(from.file){
    ad <- sapply(json, function(i) i$sampleinfo[[1]]$AD)
    names <- as.character(sapply(json, function(i) i$sampleinfo[[1]]$NAME))
    gt <- sapply(json, function(i) i$sampleinfo[[1]]$GT)
  } else {
    ad <- sapply(json$sampleinfo, function(i) i$AD)
    names <- sapply(json$sampleinfo, function(i) i$NAME)
    gt <- sapply(json$sampleinfo, function(i) i$GT)
  }
  null.idx <- sapply(ad, is.null)
  if(any(null.idx)) ad[which(null.idx)] <- '0,0'
  ad <- unlist(ad)
  vr <- VRanges(seqnames=as.character(if(from.file) sapply(json, function(i) i$chr) else json$chr), 
                ranges = IRanges(start=as.integer(if(from.file) sapply(json, function(i) i$pos) else json$pos), 
                                 end=as.integer(if(from.file) sapply(json, function(i) i$pos) else json$pos)),
                ref=as.character(if(from.file) sapply(json, function(i) i$ref) else json$ref),
                alt=as.character(if(from.file) sapply(json, function(i) i$alt) else json$alt),
                totalDepth=sapply(strsplit(ad, ','), function(i) sum(as.integer(i))),
                refDepth=as.integer(sapply(strsplit(ad, ','), function(i) i[[1]])),
                altDepth=as.integer(sapply(strsplit(ad, ','), function(i) i[[2]])),
                sampleNames = names)
  vr$GT <- as.character(gt)
  return(vr)
}