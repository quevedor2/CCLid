% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapVcf2Affy.R
\name{mapVcf2Affy}
\alias{mapVcf2Affy}
\title{mapVcf2Affy}
\usage{
mapVcf2Affy(vcfFile, snp6.dat)
}
\arguments{
\item{vcfFile}{Absolute path to VCF file}

\item{snp6.dat}{SNP6 probeset genomic position, accessible from CCLid::ccl_table}
}
\value{
List composed of two dataframes, BAF and GT
}
\description{
Takes a VCF files and reads it into VariantAnnotation package to overlap with
GRanges probeset dataset of snp6.dat.  It willc onvert the 0/0, 0/1, etc.. genotypes to 0,1,2.
Finally, it will reutnr the BAFs and normalized BAFs for each of the overlapping probesets.
}
