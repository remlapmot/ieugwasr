% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afl2.r
\name{afl2_list}
\alias{afl2_list}
\title{Retrieve a allele frequency and LD scores for pre-defined lists of variants}
\usage{
afl2_list(variantlist = "reduced", opengwas_jwt = get_opengwas_jwt(), ...)
}
\arguments{
\item{variantlist}{Choose pre-defined list. \code{"reduced"} = ~20k SNPs that are
common in all super populations (default). \code{"hapmap3"} = ~1.3 million hm3 SNPs}

\item{opengwas_jwt}{Used to authenticate protected endpoints. Login to \url{https://api.opengwas.io} to obtain a jwt. Provide the jwt string here, or store in .Renviron under the keyname OPENGWAS_JWT.}

\item{...}{Additional arguments passed to \code{api_query()}.}
}
\value{
Data frame containing ancestry specific LD scores and allele frequencies for each variant
}
\description{
Data frame includes 1000 genomes metadata including sample sizes,
allele frequency and LD score, separated by 5 super populations
(EUR = European, AFR = African, EAS = East Asian, AMR = Admixed American,
SAS = South Asian)
}
