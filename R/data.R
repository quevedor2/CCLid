#' BED file of SNP6 Probeset IDs and genomic loci in hg19
#' @docType data
#'
#' @usage data(snp6.dat)
#' @keywords datasets
#' 
#' @format A list containing 3 elements, each containing 'ProbesetID, chrom, pos, strand'
#' \describe{
#'   \item{SNP}{SNP probeset data-frame containing 933422 rows x 4 cols}
#'   \item{CN}{CN probeset data-frame containing 945615 rows x 4 cols}
#'   \item{All}{SNP and CN probeset data-frame containing 1879037 rows x 4 cols}
#'   ...
#' }
#' @source \url{http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/GenomeWideSNP_6.hg19.bed.zip}
"snp6.dat"

#' Mapping of cell IDs to filenames
#' @docType data
#'
#' @usage data(meta.df)
#' @keywords datasets
#' 
#' @format A data.frame structure
#' \describe{
#'   \item{ID}{Unique cell name}
#'   \item{GDSC}{Filename ID for GDSC data}
#'   \item{CCLE}{Filename ID for CCLE data}
#'   ...
#' }
"meta.df"