#' driftOverlap
#'
#' @param seg Seg data
#' @param ref.ds Reference dataset
#' @param alt.ds Alternate dataset
#' @importFrom methods as
#' @importFrom stats setNames
#' @importFrom utils combn
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges seqnames
#' @export
driftOverlap <- function(seg, ref.ds=NULL, alt.ds=NULL){
  if(is.null(ref.ds)) stop("Requires input of ref.ds (GDSC or CCLE)")
  if(is.null(alt.ds)) stop("Requires input of alt.ds (GDSC or CCLE)")
  
  # seg <- seg.sig[-null.idx][[1]]
  if(is.null(unlist(seg))) {
    na.mat <- matrix(rep(NA, 4), ncol=1)
    return(list(na.mat, na.mat, na.mat))
  }
  grl <- as(unlist(seg), "GRangesList")
  grl.idx <- setNames(c(.grepNA(paste0(ref.ds, "_.*", alt.ds), names(grl)),
                        .grepNA(paste0("RNA_.*", alt.ds), names(grl)),
                        .grepNA(paste0("RNA_.*", ref.ds), names(grl))),
                      c(paste0(ref.ds, "/", alt.ds),
                        paste0("RNA/", alt.ds),
                        paste0("RNA/", ref.ds)))
  
  drift.ov <- apply(combn(x=1:3, m=2), 2, function(i){
    # i <- unlist(cs[,2])
    gr1 <- if(is.na(grl.idx[i[1]]))  .blankGr() else grl[[grl.idx[i[1]]]]
    gr2 <- if(is.na(grl.idx[i[2]]))  .blankGr() else grl[[grl.idx[i[2]]]]
    grI <- intersect(gr1, gr2)
    
    wI <- sum(width(grI))
    wG1 <- sum(width(setdiff(gr1, grI)))
    wG2 <- sum(width(setdiff(gr2, grI)))
    wEmpty <- 0
    
    gr1.chr <- as.character(seqnames(gr1))
    gr2.chr <- as.character(seqnames(gr2))
    if(any(gr1.chr == 'chrZ') & any(gr2.chr == 'chrZ')){
      wEmpty <- 1
      wI <- 0
    }
    
    w.frac <- as.matrix(round(c(wG1, wI, wG2, wEmpty) / sum(c(wG1, wI, wG2, wEmpty)),2))
    rownames(w.frac) <- c(names(grl.idx)[i[1]], 'intersect', names(grl.idx)[i[2]], 'no_drift')
    return(as.data.frame(w.frac))
  })
  return(drift.ov)
}
