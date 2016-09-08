#' Make a population
#' 
#' @description 
#' Takes haploid genome information and a crossing block and produce a list of 
#' families. The function also requires the genome, the size of the families 
#' and the mutation rates.
#' 
#' @param genome The hypred genome.
#' @param parental.haploids The complete 2n x m matrix of parent haploid
#' genotype information.
#' @param N The number of F_1 individuals to generate.
#' @param generations The number of selfing generations of single-seed descent. 
#' For instance \code{generations = 5} would results in an F_6 population. 
#' For each generation, each of the F_1 individuals results in a single new 
#' individual that has undergone an additional generation of inbreeding
#' @param pop.type The type of population to generate. Can be \code{"inbred"} 
#' or \code{"random"} (default).
#' @param cycle.number The integer number of the breeding cycle.
#' @param mutate A logical as to generate random mutations when recombining.
#' @param mutation.rate.snp The per-base mutation rate of the SNPs.
#' @param mutation.rate.qtl The per-base mutation rate of the QTL.
#' 
#' @import dplyr
#' @import stringr
#' 
#' @export
#' 
make.population <- function(genome, parental.haploids, crossing.block, N, 
                            cycle.number, generations, pop.type = "random",
                            mutate = TRUE, mutation.rate.snp, mutation.rate.qtl) {
  
  # Make sure the crossing block is a matrix
  crossing.block <- as.matrix(crossing.block)
  # Pull out the loci names
  loci.names <- colnames(parental.haploids)
  
  # The number of crosses
  n.crosses = nrow(crossing.block)
  
  # Create empty vector
  fam.list <- vector("list", n.crosses)
  
  # Iterate over the number of crosses
  for (i in seq(n.crosses)) {
    
    # Extract the cross
    cross <- crossing.block[i,]
    
    # Pull out the gametes
    parent1.haploid <- parental.haploids %>%
      subset.matrix(subset = str_detect(string = row.names(.), pattern = cross[1]))
    parent2.haploid <- parental.haploids %>%
      subset.matrix(subset = str_detect(string = row.names(.), pattern = cross[2]))
    
    # Make a population from each cross
    pop <- make.family(genome = genome, 
                       parent1.genome = parent1.haploid, 
                       parent2.genome = parent2.haploid, 
                       N = N, 
                       generations = generations, 
                       pop.type = pop.type,
                       family.number = i,
                       cycle.number = cycle.number,
                       mutate = mutate,
                       mutation.rate.snp = mutation.rate.snp, 
                       mutation.rate.qtl = mutation.rate.qtl)
    
    # Add marker names back in
    colnames(pop) <- loci.names
    
    # Add the population to the list
    fam.list[[i]] <- pop
    
  } # Close the per-cross loop
  
  # Rename the list entries
  names(fam.list) <- apply(X = crossing.block, MARGIN = 1, FUN = function(cross) return(paste(cross[1],cross[2], sep = ".")))
  
  # Return the list
  return(fam.list)
  
} # Close the function