#' Measure expected heterozygosity
#' 
#' @description 
#' Calculated expected heterozygosity as the mean expected heterozygosity
#' across all m loci
#' 
#' @param genome The hypred genome
#' @param haploid.genos The 2n x m matrix of haploid genotypes for n entries
#' and m SNPs and QTL.
#' 
#' @import dplyr
#' 
#' @export
#' 
measure.expected.het <- function(genome, haploid.genos) {
  
  # Deal with input
  if(!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- do.call("rbind", haploid.genos)
  }
  
  # Calculate the allele frequency of the loci
  allele.freq <- measure.af(genome = genome, haploid.genos = haploid.mat)$loci
  # Combine with minor allele frequency
  freqs <- data.frame(allele1 = allele.freq, 
                      allele0 = 1 - allele.freq)
  
  # Calculate the sum of the squared allele frequencies at each locus
  summation <- (freqs^2) %>% rowSums()
  
  # Average the summation and subtract that from 1
  1 - mean(summation)
  
} # Close the function