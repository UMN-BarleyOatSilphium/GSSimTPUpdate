#' Find positions of loci
#' 
#' @description
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