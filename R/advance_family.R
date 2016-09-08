#' Advance a family
#' 
#' 
#' 
#' 
advance.family <- function(genome, 
                           pop.mat, 
                           pop.type = "random",
                           cycle.number,
                           starting.generation, 
                           generations, 
                           mutation.rate.snp, 
                           mutation.rate.qtl) {
  
  # Deal with inputs
  pop.mat <- as.matrix(pop.mat)
  # Gather marker names
  markers <- colnames(pop.mat)
  
  # Rename the population
  pop.i <- pop.mat
  
  # Split line information on the inbreeding generation
  line.names <- row.names(pop.i)
  line.name.length <- unique(sapply(X = line.names, FUN = nchar))
  
  prefix <- paste("C", cycle.number, sep = "")
  suffix <- substring(text = line.names, first = (line.name.length - 9))
  
  # Find the number of gametes
  n.gamete <- nrow(pop.i)
  
  # If no advancement is requested, don't advance
  if (generations > 0) {
    # Iterate over generations
    for (g in 1:generations) {
      
      # This applies the function over individuals
      # First we generate a set of recombinant gametes
      recomb.pop <- lapply(X = seq(1, n.gamete, 2), FUN = function(i) {
        gamete.index1 <- i
        gamete.index2 <- i + 1
        
        # Simulate the gametes from the i-th individual
        gamete1 <- 
          hypredRecombine(genome,
                          genomeA = pop.i[gamete.index1,],
                          genomeB = pop.i[gamete.index2,],
                          mutate = T,
                          mutation.rate.snp = mutation.rate.snp,
                          mutation.rate.qtl = mutation.rate.qtl,
                          block = FALSE)
        
        gamete2 <-
          hypredRecombine(genome,
                          genomeA = pop.i[gamete.index1,],
                          genomeB = pop.i[gamete.index2,],
                          mutate = T,
                          mutation.rate.snp = mutation.rate.snp,
                          mutation.rate.qtl = mutation.rate.qtl,
                          block = FALSE)
        
        return( list(gamete1, gamete2) )
      })
      
      # Collapse the list into the correct matrix
      pop.i <- do.call("rbind", lapply(X = recomb.pop, FUN = function(x) do.call("rbind", x)))
      
      # If randomly mating, permutate
      if (pop.type == "random") {
        pop.i <- pop.i[sample(1:n.gamete),]
      }
      
    } # Close the generation for loop
  } # Close the if statement
  
  # Assign line names
  row.names(pop.i) <- paste(prefix, "_", (starting.generation + generations), suffix, sep = "")
  # Assign marker names
  colnames(pop.i) <- markers
  
  # Return the population matrix
  return(pop.i)
} # Close the function