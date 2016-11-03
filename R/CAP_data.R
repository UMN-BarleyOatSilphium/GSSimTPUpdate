#' CAP genotypic data in rrBLUP encoding
#' 
#' @description
#' A \code{data.frame} of genotype data in rrBLUP encoding. The column
#' names are metadata variables and line names. This genotype data is unfiltered. 
#' 
#' @format A data frame with 2832 rows and 768 columns:
#' \describe{
#'    \item{rs}{The SNP marker name}
#'    \item{alleles}{The SNP marker alleles}
#'    \item{chrom}{The chromosome in which the marker resides}
#'    \item{pos}{The genetic position (in thousandths of cM)}
#' }
#' 
#' @source \url{https://triticeaetoolbox.org/barley/}
"CAP.hmp"


#' CAP genotypic data in haploid matrix format
#' 
#' @description 
#' A \code{matrix} of genotype data in haploid format. The column
#' names are the names of the SNP markers and the row names are the names 
#' of the lines. Each line has two names with the suffix "1" or "2" to denote
#' the haplotype number. These haploids were created using the 
#' "CAP_data_preparation.R" script.
#' 
#' @format A data frame with 1528 rows and 1590 columns
#' 
#' @source \url{https://triticeaetoolbox.org/barley/}
"CAP.haploids"

#' The names of the CAP lines
#' 
#' @description 
#' A \code{character} of CAP line names for which genotypic data was 
#' used in the simulations.
#' 
#' @format A character with 764 elements.
#' 
#' @source \url{https://triticeaetoolbox.org/barley/}
"CAP.lines"

#' Marker information
#' 
#' @description 
#' A \code{data.frame} of marker metadata.
#' 
#' @format A data frame with 1590 rows and 4 columns
#' \describe{
#'    \item{rs}{The SNP marker name}
#'    \item{alleles}{The SNP marker alleles}
#'    \item{chrom}{The chromosome in which the marker resides}
#'    \item{pos}{The genetic position (in Morgans)}
#' }
#' 
#' @source \url{https://triticeaetoolbox.org/barley/}
"CAP.markers"
