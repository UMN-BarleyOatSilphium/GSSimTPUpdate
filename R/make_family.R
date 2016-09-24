#' Make a family
#' 
#' @param genome The list of hypred genomes.
#' @param parent1.genome A 2 x m matrix of the parent 1 haploid genomes.
#' @param parent2.genome See \code{parent1.genome}
#' @param N The number of F_1 individuals to generate.
#' @param generations The number of selfing generations of single-seed descent. 
#' For instance \code{generations = 5} would results in an F_6 population. 
#' For each generation, each of the F_1 individuals results in a single new 
#' individual that has undergone an additional generation of inbreeding
#' @param pop.type The type of population to generate. Can be \code{"inbred"} 
#' or \code{"random"} (default).
#' @param cycle.number The integer number of the breeding cycle.
#' @param family.number The integer number of the family.
#' @param mutate A logical as to generate random mutations when recombining.
#' @param mutation.rate.snp The per-base mutation rate of the SNPs.
#' @param mutation.rate.qtl The per-base mutation rate of the QTL.
#' 
#' 
#' @import dplyr
#' 
#' @export
#' 
#' 
make.family <- function(genome, parent1.genome, parent2.genome, N, generations, 
                        pop.type = "random", cycle.number, family.number, mutate = TRUE,
                        mutation.rate.snp, mutation.rate.qtl) {

  # Determine the number of gametes to produce
  n.gamete = N * 2
  
  
  F_1 <- replicate(n = n.gamete, expr = {
    
    rbind(
      # Gamete 1
      recombine(genome = genome,
                haploid.genomeA = parent1.genome[1,],
                haploid.genomeB = parent1.genome[2,],
                mutate = mutate,
                mutation.rate.snp = mutation.rate.snp,
                mutation.rate.qtl = mutation.rate.qtl),
      # Gamete 2
      recombine(genome = genome,
                haploid.genomeA = parent2.genome[1,],
                haploid.genomeB = parent2.genome[2,],
                mutate = mutate,
                mutation.rate.snp = mutation.rate.snp,
                mutation.rate.qtl = mutation.rate.qtl) ) %>%
      list() }) %>%
    # Transpose
    do.call("rbind", .)
  
  # Iterate changes to the population over a number of generations
  # Rename the population
  pop.i <- F_1
  
  # If no advancement is requested, don't advance
  if (generations > 0) {
    
    # Iterate over generations
    for (g in seq(generations)) {
      
      # This applies the function over individuals
      # First we generate a set of recombinant gametes
      pop.i <- lapply(X = seq(1, n.gamete, 2), FUN = function(i) {
        gamete.index1 <- i
        gamete.index2 <- i + 1
        
        # Simulate the gametes from the i-th individual
        replicate(n = 2, expr = {
          
          recombine(genome = genome,
                    haploid.genomeA = pop.i[gamete.index1,],
                    haploid.genomeB = pop.i[gamete.index2,],
                    mutate = mutate,
                    mutation.rate.snp = mutation.rate.snp,
                    mutation.rate.qtl = mutation.rate.qtl) }) %>%
          # Transpose
          t() }) %>%
        
        # Bind the rows
        do.call("rbind", .)
      
      # If randomly mating, permutate
      if (pop.type == "random") {
        pop.i <- pop.i[sample(1:n.gamete),]
      }
      
    } # Close the generation for loop
    
  } # Close the if statement
  
  # Name the individuals
  row.names(pop.i) <- paste("C", cycle.number, "_", (1+generations), 
                            formatC(family.number, width = 3, format = "d", flag = "0"), 
                            "-", formatC(rep(c(1:N), each = 2), width = 4, format = "d", flag = "0"),
                            ".", rep(c(1:2), 2), sep = "")
  
  # Return the population matrix
  return(pop.i)
  
} # Close the function