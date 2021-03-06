% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demo_data.R
\name{combineSamples}
\alias{combineSamples}
\title{combineSamples}
\usage{
combineSamples(data.type, sample.mat, prop)
}
\arguments{
\item{data.type}{'BAF' or 'geno'}

\item{sample.mat}{Input matrix of BAF or geno values by samples}

\item{prop}{Proportion to mix the 2+ samples by}
}
\value{
A weighted mean vector of BAF or a probability sampled genotype
}
\description{
Combines the BAF or Geno or two or samples given a set proportion
}
