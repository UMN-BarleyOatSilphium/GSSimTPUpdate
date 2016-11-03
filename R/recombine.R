#' Simulate recombination
#' 
#' @description 
#' Uses the \code{\link[hypred]{hypredRecombine}} function as a base, but will 
#' apply that function over the chromosomes in the genome list
#' 
#' @param genome The list of hypred genomes.
#' @param haploid.genomeA The first haploid genome.
#' @param haploid.genomeB The second haploid genome.
#' @param mutate A logical as to generate random mutations when recombining.
#' @param mutation.rate.snp The per-base mutation rate of the SNPs.
#' @param mutation.rate.qtl The per-base mutation rate of the QTL.
#' 
#' @import hypred
#' 
#' @export
#' 
recombine <- function(genome, haploid.genomeA, haploid.genomeB, mutate = TRUE,
                      mutation.rate.snp, mutation.rate.qtl) {
  
  # Pull out the number of chromosomes
  n.chr <- length(genome)
  # Pull out the number of loci per chromsosome
  loci.per.chr <- sapply(X = genome, FUN = function(chr) length(chr@pos.snp))
  
  # Split each haploid genome into chromosomes
  haploid.genomeA.split <- split(haploid.genomeA, rep(seq(n.chr), loci.per.chr))
  haploid.genomeB.split <- split(haploid.genomeB, rep(seq(n.chr), loci.per.chr))
  
  # Apply a function over each chromosome
  gamete.list <- mapply(genome, haploid.genomeA.split, haploid.genomeB.split,
                        FUN = function(chr, hgA, hgB) {
    
    hypredRecombine(chr,
                    genomeA = hgA,
                    genomeB = hgB,
                    mutate = mutate,
                    mutation.rate.snp = mutation.rate.snp,
                    mutation.rate.qtl = mutation.rate.qtl,
                    block = FALSE) })
  
  # Concatenate the gamete
  return(do.call("c", gamete.list))
  
} # Close the function