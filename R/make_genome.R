#' Make a genome
#' 
#' @param n.chr The number of chromosomes in the genome
#' @param chr.len A vector of length \code{n.chr} with elements corresponding to
#' the Morgan length of each chromosome.
#' @param n.chr.snps A vector of length \code{n.chr} with elements corresponding
#' to the number of SNP loci on each chromosome.
#' @param genetic.map A list of length \code{n.chr} where each element is a vector
#' of Morgan positions corresponding to the SNP loci.
#' 
#' @return A \code{list} of \code{hypred} genomes.
#' 
#' @import hypred
#' 
#' @export
#' 
make.genome <- function(n.chr, chr.len, n.chr.snps, genetic.map) {
  
  # Apply a function over the number of chromosomes
  genome.list <- lapply(X = 1:n.chr, FUN = function(i)
    # Create a new base genome
    hypredGenome(num.chr = 1, len.chr = chr.len[i], num.snp.chr = n.chr.snps[i]) )
  
  # Rename the genomes
  names(genome.list) <- paste("chr", 1:n.chr, sep = "")
  
  # Apply a function over the genome list to add a new genetic map
  genome.list <- lapply(X = 1:n.chr, FUN = function(i) {
    # Extract the genome
    genome <- genome.list[[i]]
    # Add the new map to the genome
    hypredNewMap(genome, new.map = genetic.map[[i]]) })
  
  # Return the genome
  return(genome.list)
  
} # Close the function