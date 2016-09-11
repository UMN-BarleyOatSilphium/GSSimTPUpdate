#' Hidden functions
#' 
#' @description 
#' Attempts to find the inverse of a matrix. If unsuccessful, find the
#' generalized inverse.
#' 
#' @param x A matrix to be inverted.
#' @param silent Logical whether to print the method of inversion used.
#' 
#' @import MASS
#' 
invert.mat <- function(x, silent = FALSE) {
  
  if (try(expr = solve(x), silent = T) %>% class == "try-error") {
    x.inv <- ginv(x)
    method <- "Generalized"
  } else {
    # Otherwise use solve
    x.inv <- solve(x)
    method <- "Inverse"
  }
  
  if (!silent)  cat("\nMethod: ", method)
  return(x.inv)
  
} # Close the function

#' 
#' Determine if a locus is polymorphic
#' 
#' @import dplyr
#' 
is.polymorphic <- function(x) {
  mean(x) %>%
    abs() != 1
} # Close the function

#' 
#' Make a matrix of contrasts
#' 
#' @description 
#' Creates a matrix of contrasts between the individuals not in the training 
#' population (i.e. the candidates) and the mean of the whole population.
#' 
#' @param unphenotyped.index The index of the unphenotyped entries in the 
#' relationship matrix (A)
#' @param n.total The size of the whole population
#' 
#' @references 
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P., 
#' Moreau, L. (2012). Maximizing the Reliability of Genomic Selection by 
#' Optimizing the Calibration Set of Reference Individuals: Comparison of 
#' Methods in Two Diverse Groups of Maize Inbreds (Zea mays L.). Genetics, 
#' 192(2), 715–728. http://doi.org/10.1534/genetics.112.141473
#' 
make.contrast <- function(unphenotyped.index, n.total) {
  
  # Number of unphenotyped individuals
  n.unphenotyped <- length(unphenotyped.index)
  
  # The positive value of the contrast
  pos.contrast = 1 - (1 / n.total)
  # The negative value of the contrast
  neg.contrast = -1 / n.total
  
  # Make the total constrast matrix
  contrast.mat <- matrix(data = neg.contrast, nrow = n.total, ncol = n.unphenotyped)
  
  # Fill in the positive values
  coord <- cbind(unphenotyped.index,
                 seq(n.unphenotyped))
  
  contrast.mat[coord] <- pos.contrast
  
  return(contrast.mat)
}

#' 
#' Finds the whole genome positions of QTL and markers
#' 
#' @param genome The list of hypred genomes.
#' @param genos The n x m matrix of genotypes for n entries across m loci.
#' 
#' @import dplyr
#' 
find.pos <- function(genome, genos) {
  
  # Return all the loci (this is easy - just the index of the number of columns)
  pos.loci <- genos %>%
    ncol() %>%
    seq()
  
  # Number of chromosomes
  n.chr <- length(genome)
  
  # Running count of the number of loci up to the nth chromosome
  n.loci.total <- 0
  
  # Empty lists of qtl and snp positions
  pos.qtl <- list()
  
  # Iterate over chromosomes to find the position of markers and qtl
  for (i in seq(n.chr)) {
    
    # Extract the position of the qtl and add the loci total
    pos.qtl[[i]] <- genome[[i]]@pos.add.qtl$ID + n.loci.total
    
    # Total number of loci up to the chromosome
    n.loci.total = n.loci.total + genome[[i]]@num.snp.chr
    
  }
  
  # Unlist the qtl position list
  pos.qtl <- pos.qtl %>%
    unlist()
  
  # Position of markers
  pos.snp <- setdiff(pos.loci, pos.qtl)
  
  # Create list
  out.list <- list(pos.loci = pos.loci,
                   pos.qtl = pos.qtl,
                   pos.snp = pos.snp)
  
  # Return
  return(out.list)
  
} # Close the function

#' 
#' Round with a limit
#'
#' @description 
#' A convenience function to round a number down or up if it is outside a 
#' lower or upper limit
#' 
round.limit <- function(x, limit, option) {
  if (option == "upper") {
    if (x > limit) {
      return(limit)
    } else {return(x) }
  }
  if (option == "lower") {
    if (x < limit) {
      return(limit)
    } else {return(x) }
  }
} # Close

#' 
#' M matrix
#' 
#' @description 
#' Creates an M matrix, or the orthagonal projector
#' 
#' @param n The number of rows of matrix X
#' 
#' @references 
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P., 
#' Moreau, L. (2012). Maximizing the Reliability of Genomic Selection by 
#' Optimizing the Calibration Set of Reference Individuals: Comparison of 
#' Methods in Two Diverse Groups of Maize Inbreds (Zea mays L.). Genetics, 
#' 192(2), 715–728. http://doi.org/10.1534/genetics.112.141473
#' 
design.M <- function(n) {
  
  # Identity matrix of size n
  I <- diag(n) 
  # X matrix
  X <- matrix(1, n, 1)
  # Create M and return
  I - tcrossprod(X %*% solve(crossprod(X, X)), X)
  
} # Close the function