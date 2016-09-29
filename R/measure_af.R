#' Calculate allele frequencies
#' 
#' @description 
#' Calculates the frequency of alleles at SNP markers or QTL. Outputs the 
#' minor allele frequency if requested.
#' 
#' @param genome The hypred genome.
#' @param genos The 2n x m matrix of haploid genotypes for n entries and m SNPs
#' and QTL.
#' 
#' @return 
#' A vector of allele frequencies for m loci.
#' 
#' @export 
#' 
measure.af <- function(genome, haploid.genos) {
  
  # Deal with input
  if(!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- do.call("rbind", haploid.genos)
  }
  
  # Find the position of QTL and SNP markers
  pos.snp <- GSsim.TPUpdate:::find.pos(genome)
  
  # The mean of each column gives us the frequency of the 1 allele at a site.
  allele.freq <- colMeans(haploid.mat)
  
  # Separate by SNP and QTL
  freq.by.type <- lapply(X = pos.snp, FUN = function(cat) 
    allele.freq[cat] )
  
  # Rename
  names(freq.by.type) <- c("loci", "qtl", "snp")
  
  return(freq.by.type)
  
  
} # Close the function