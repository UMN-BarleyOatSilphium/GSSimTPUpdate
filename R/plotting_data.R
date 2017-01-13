#' Results
#' 
#' @description
#' A \code{list} of simulation results. Each list element is a different variable
#' measured over the course of the simulations.
#' 
#' @format A list with 9 measured variables:
#' \describe{
#'    \item{accuracy}{The prediction accuracy}
#'    \item{genvar}{Genetic variance}
#'    \item{genval}{Genotypic values}
#'    \item{qtl_marker_LD}{Different measures of linkage disequilibrium}
#'    \item{rel}{Average genomic relationship}
#'    \item{inbreeding}{Average inbreeding of the selection candidates}
#'    \item{fixedqtl}{The proportion of QTL fixed for an allele}
#'    \item{exphet}{Average expected heterozygosity measured on the parents}
#'    \item{fixedmar}{The proportion of markers fixed for an allele}
#' }
#' 
"plotting_data"