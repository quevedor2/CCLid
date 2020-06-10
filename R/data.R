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

#' Mapping of Omni 2.5M probesets to Affy SNP6.0
#' @docType data
#'
#' @usage data(affy.omni)
#' @keywords datasets
#' 
#' @format An annotated data.frame structure
#' \describe{
#'   \item{mergeid}{Unique snp identifier}
#'   \item{Probe_Set_ID}{Affy6 Probeset ID}
#'   \item{dbSNP_RS_ID}{dbSNP RS ID}
#'   \item{Chromosome}{Affy6 Chr}
#'   \item{Physical_Position}{Affy6 genomic position}
#'   \item{Strand}{Affy6 Strand}
#'   \item{Allele_A}{Affy6 Allele A}
#'   \item{Allele_B}{Affy6 Allele B}
#'   \item{Name}{Omni2.5 ID}
#'   \item{SNP}{Omni2.5 SNPs in X/Y format}
#'   \item{Chr}{Omni2.5 Chromosome}
#'   \item{MapInfo}{Omni2.5 hg19 genomic position}
#'   \item{RefStrand}{Omni2.5 Strand}
#'   \item{SNP_A}{Omni2.5 Allele A for RefStrand}
#'   \item{SNP_B}{Omni2.5 Allele B for RefStrand}
#'   \item{scSNP_A}{Omni2.5 strand-corrected Allele A}
#'   \item{scSNP_B}{Omni2.5 strand-corrected Allele B}
#'   \item{flip}{Value to flip Omni2.5 alleles if strands didnt match}
#'   ...
#' }
"affy.omni"

#' Chromosomal instability (CIN) 70 genes from https://doi.org/10.1038/ng1861
#' @docType data
#'
#' @usage data(cin70)
#' @keywords datasets
#' 
#' @format A data.frame structure
#' \describe{
#'   \item{Gene_old}{Published gene names from the paper, some are outdated}
#'   \item{Gene}{HUGO symbol ID updates for the old published genes}
#'   \item{core}{Whether the genes are core CIN genes (1) or not (0)}
#'   \item{ENS}{Ensemble ENS ID number}
#'   ...
#' }
"cin70"

#' Metadata for gCSI dataset
#' @docType data
#'
#' @usage data(gne.meta)
#' @keywords datasets
#' 
#' @format A data.frame structure
#' \describe{
#'   ...
#' }
"gne.meta"

#' Metadata mapping together all pharmacogenomic datasets and their RNA IDs
#' @docType data
#'
#' @usage data(rna.meta.df)
#' @keywords datasets
#' 
#' @format A data.frame structure
#' \describe{
#'   ...
#' }
"rna.meta.df"

#' Metadata mapping together all pharmacogenomic datasets and their RNA IDs
#' @docType data
#'
#' @usage data(meta.df)
#' @keywords datasets
#' 
#' @format A data.frame structure
#' \describe{
#'   ...
#' }
"meta.df"
