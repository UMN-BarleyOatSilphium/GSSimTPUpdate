#' Advance a population
#' 
#' @description
#' Takes a list of populations from a cross and advance each population either 
#' through random mating or through inbreeding.
#' 
#' 
#' 
advance.population <- function(genome, 
                               list.of.families, 
                               pop.type = "random", 
                               generations,
                               cycle.number,
                               starting.generation,
                               mutation.rate.snp, 
                               mutation.rate.qtl) {
  
  # Apply a function over the list
  advanced.family <- lapply(X = list.of.families, FUN = function(family) {
    advanced.family <- advance.family(genome = genome, 
                                      pop.mat = family, 
                                      pop.type = pop.type, 
                                      generations = generations, 
                                      cycle.number = cycle.number, 
                                      starting.generation = starting.generation,
                                      mutation.rate.snp = mutation.rate.snp,
                                      mutation.rate.qtl = mutation.rate.qtl)
    return(advanced.family)
  })
  
  return(advanced.family)
} # Close the function
