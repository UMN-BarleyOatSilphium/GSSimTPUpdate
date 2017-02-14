#' Results
#' 
#' @description
#' A \code{list} of simulation results. Each list element is a different variable
#' measured over the course of the simulations.
#' 
#' @format A list with 11 \code{data.frame} objects with a different measured variable:
#' \describe{
#'    \item{df.acc}{The prediction accuracy}
#'    \item{df.genvar}{Genetic variance}
#'    \item{df.genval}{Genotypic values}
#'    \item{df.resp}{Response to selection}
#'    \item{df.tpmeanmax}{QTL-marker LD in the training population, measured as the average max LD of a QTL with any marker across the genome}
#'    \item{df.scmeanmax}{QTL-marker LD in the selection candidates, measured as the average max LD of a QTL with any marker across the genome}
#'    \item{df.pers}{Persistence of LD phase}
#'    \item{df.rel}{Average genomic relationship}
#'    \item{df.inbred}{Average inbreeding of the selection candidates}
#'    \item{df.rateinbred}{Average rate of inbreeding in the selection candidates}
#'    \item{df.fixedqtl}{The proportion of QTL fixed for an allele}
#' }
#' 
"plotting_data"