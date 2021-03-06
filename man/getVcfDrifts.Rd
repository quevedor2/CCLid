% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paper_functions.R
\name{getVcfDrifts}
\alias{getVcfDrifts}
\title{getVcfDrifts}
\usage{
getVcfDrifts(
  vcfFile,
  ref.dat,
  meta.df,
  min.depth = 5,
  centering = "extreme",
  dataset = "GDSC",
  snp6.dat
)
}
\arguments{
\item{vcfFile}{path to vcf file}

\item{ref.dat}{Refrence data containing Reference matrix and variance}

\item{meta.df}{Meta file linking RNA files to cell names}

\item{min.depth}{Minimum depth to consider for SNPs}

\item{centering}{Centering of data method to be passed into bafDrift() function}

\item{dataset}{Either 'CCLE', 'GDSC', ro 'GNE'}

\item{snp6.dat}{SNP6 probeset genomic position, accessible from CCLid::ccl_table}
}
\value{
A list containing the fraction of genome drifted,
as well the significantly drifted regions CNAo
}
\description{
Check the drift of a given VCF file and all matching
cell line names based on the meta-data matched IDs
}
