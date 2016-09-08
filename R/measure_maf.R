#' Calculate the minor allele frequency across a genotype matrix
#' 
#' @description 
#' Calculates the minor allele frequency. Also outputs the site-frequency
#' spectrum, if desired.
#' 
#' @param genos The n x m matrix of genotypes for n entries across m loci.
#' @param plot.sfs A logical whether to plot the site-frequency spectrum.
#' 
#' @return 
#' A vector of minor allele frequencies for m loci.
#' 
#' @export 
#' 
measure.maf <- function(genos, plot.sfs = FALSE) {
  
  # Deal with input
  genos <- as.matrix(genos)
  
  # Add 1 to get the count of the 1 allele at each locus per entry
  genos1 <- genos + 1
  
  # Taking the mean across a SNP gives the frequency of the 1 allele as a function
  ## of diploid genotypes. We need to divide by two to get the frequency as a function
  ## of haploid genotypes
  freq <- colMeans(genos1) / 2
  maf <- sapply(X = freq, FUN = function(f) min(f, 1-f))
  
  
  # If requested, plot the site-frequency spectrum
  if(plot.sfs) {
    
    # Plot
    hist(x = maf,
         main = "Site Frequency Spectrum",
         xlab = "Minor Allele Frequency",
         ylab = "Count",
         xlim = c(0,0.5))
  }
  
  # Return the data
  return(maf)
  
} # Close the function