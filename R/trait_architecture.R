#' Define the trait architecture
#' 
#' @description 
#' Adds QTL to a genome. SNPs are assigned (randomly or by the user) to become
#' QTL. Effects are then assigned to the QTL (from various distributions). All 
#' QTL have additive effects, and some, all, or none can have dominance effects.
#' 
#' @param genome The list of hypred genomes.
#' @param n.QTL The number of SNPs to become QTL.
#' @param qtl.index A list of SNP indices per chromosome to become QTL. If NULL, 
#' indices are randomly sampled.
#' @param qtl.dom.index A list of QTL indices per chromosome to have dominance 
#' effects.
#' @param qtl.perf.index A list of QTL indices per chromosome to be perfect
#' markers (i.e. the observed SNP is the QTL).
#' @param qtl.add.eff A vector of length \code{n.QTL} of additive effects of QTL.
#' The vector is interpreted as such: the ith element in the vector will be the
#' additive effect assigned to the ith QTL. Can also be \code{"normal"} to draw 
#' effects from a standard normal distribution or \code{"geometric"} to draw
#' effects from a geometric series.
#' @param qtl.dom.eff A vector of length \code{n.QTL} of dominance effects of QTL.
#' The vector is interpreted as such: the ith element in the vector will be the
#' dominance effect assigned to the ith QTL.
#' 
#' @return 
#' A list of hypred genomes with assigned trait architecture.
#' 
#' @import hypred
#' 
#' @export
#' 
trait.architecture <- function(genome, n.QTL, qtl.index = NULL, 
                               qtl.dom.index = NULL, qtl.perf.index = NULL, 
                               qtl.add.eff = "normal", qtl.dom.eff = NULL) {
  
  # Deal with input
  # Find the number of chromosomes
  n.chr <- length(genome)
  
  # If the qtl.ids are NULL, randomly sample the loci
  if (is.null(qtl.index)) {
    
    # Create a list to determine the number of qtl for each chromosome. This will
    ## essentially make the distribution of qtl even across the chromosome indices
    qtl.per.chr <- sapply(X = split(x = 1:n.QTL, f = cut(1:n.QTL, breaks = n.chr)), FUN = length)
    
    # Create an empty list
    qtl.index.per.chr <- list()
    
    # Loop over the chromosomes in the genome
    for (p in 1:n.chr) {
      
      # Pull out the number of QTL designated for the chromosome
      n.qtl.p <- qtl.per.chr[p]
      
      # Pull out the number of loci on the chromsome and create an index
      loci.index.p <- seq(1, slot(genome[[p]], "num.snp.chr"))
      
      # Sample the index and add to the list
      qtl.index.per.chr[[p]] <- sort(sample(x = loci.index.p, size = n.qtl.p))
    }
    
    # If the list is specified, check it
  } else {
    
    if (length(qtl.index) != n.chr) stop("The length of the qtl.index list is not the same as the number of chromosomes.")
    # Rename it
    qtl.index.per.chr <- qtl.index
    # Figure out the QTL per chromosomes
    qtl.per.chr <- sapply(X = qtl.index.per.chr, FUN = length)
  }
  
  # QTL dominance IDs
  if (is.null(qtl.dom.index)) {
    # Create a list of NULLs
    qtl.dom.index.per.chr <- sapply(X = 1:n.chr, FUN = function(i) NULL)
    
  } else { # Otherwise check the list
    if (length(qtl.dom.index) != n.chr) stop("The length of the qtl.dom.index list is not the same as the number of chromsomes.")
    
    # Rename it
    qtl.dom.index.per.chr <- qtl.dom.index
  }
  
  # Perfect QTL IDs
  if (is.null(qtl.perf.index)) {
    # Create a list of NULLs
    qtl.perf.index.per.chr <- sapply(X = 1:n.chr, FUN = function(i) NULL)
    
  } else { # Otherwise check the list
    if (length(qtl.perf.index) != n.chr) stop("The length of the qtl.perf.index list is not the same as the number of chromosomes.")
    
    # Rename it
    qtl.perf.index.per.chr <- qtl.perf.index
  }
  
  # Assign qtl effects
  if (is.character(qtl.add.eff)) {
    
    # Draw from standard normal
    if (qtl.add.eff == "normal") {
      qtl.add.eff <- abs(rnorm(n = n.QTL, mean = 0, sd = 1))
      # Randomly assign negative values to the qtl effects - this corresponds to the value of the 1 allele
      qtl.add.eff <- qtl.add.eff * sample(c(1,-1), n.QTL, replace = T)
      
    } else { # Draw from geometric series
      
      if (qtl.add.eff == "geometric") {
        a = (n.QTL - 1) / (n.QTL + 1)
        qtl.add.eff <- sample(a^(1:n.QTL))
        # Randomly assign negative values to the qtl effects - this corresponds to the value of the 1 allele
        qtl.add.eff <- qtl.add.eff * sample(c(1,-1), n.QTL, replace = T)
        
        # Break up the effects into chromosomes with the same number of elements
        ## as QTL on those chromosomes
        qtl.add.eff.per.chr <- split(qtl.add.eff, rep(1:n.chr, qtl.per.chr))
        
      } else {
        # Otherwise, make sure the length of the additive effects is the same as the number of chr
        if (length(qtl.add.eff) != n.chr) stop("The length of the QTL effects is not the same as the number of chromosomes")
      }}}
  
  # Dominance effects
  if (is.null(qtl.dom.eff)) {
    # Create a list of NULLs
    qtl.dom.eff.per.chr <- sapply(X = 1:n.chr, FUN = function(i) NULL)
    
  } else { # Otherwise check it
    if (length(qtl.dom.eff) != n.chr) stop("The length of the qtl.dom.eff list is not the same as the number of chromosomes.")
    
    # Make sure each element in the qtl.dom.eff list is the same as the qtl.dom.index.per.chr list
    dom.elements.same <- sapply(X = 1:n.chr, FUN = function(i) length(qtl.dom.eff[[i]]) == length(qtl.dom.index.per.chr[[i]]) )
    
    if (!all(dom.elements.same)) stop("One or more of the elements in the qtl.dom.eff list are not the same length as the corresponding qtl.dom.index list.")
  }
  
  
  
  # Apply a function over each chromosome in the genome
  genome <- sapply(X = 1:n.chr, FUN = function(i) {
    
    # Add to the ith position in the genome list a revised chromosome genome with
    ## the qtl positions and effects
    genome[[i]] <- hypredNewQTL(genome[[i]],
                                new.id.add = qtl.index.per.chr[[i]],
                                new.id.dom = qtl.dom.index.per.chr[[i]],
                                new.id.per.mar = qtl.perf.index.per.chr[[i]],
                                new.eff.add = qtl.add.eff.per.chr[[i]],
                                new.eff.dom = qtl.dom.eff.per.chr[[i]] )
  })
  
  # Return the new genome
  return(genome)
  
} # Close the function
