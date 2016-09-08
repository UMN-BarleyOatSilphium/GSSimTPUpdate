#' Measure expected heterozygosity
#' 
#' @description 
#' Calculated expected heterozygosity as the mean expected heterozygosity
#' across all m loci
#' 
#' @param genos The n x m matrix of genotypes for n entries across m loci.
#' 
#' @import dplyr
#' 
#' @export
#' 
measure.expected.het <- function(genos) {
  
  # Calculate the minor allele frequency of the loci
  maf <- measure.maf(genos)
  
  # Create a matrix of the frequency of the two alleles
  allele.freqs <- cbind(maf, 1-maf)
  
  # Iterate over loci
  summation <- apply(X = allele.freqs, MARGIN = 1, FUN = function(locus)
    # Square the frequencies and sum
    locus^2 %>%
      sum() )
  
  # Average the summation and subtract that from 1
  1 - mean(summation)
  
} # Close the function