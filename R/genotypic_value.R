#' Calculate the genotypic value of individuals
#' 
#' @description 
#' Calculates the genotypic value of a set of lines based on their QTL alleles 
#' and the effect of those alleles.
#' 
#' @param genome The hypred genome.
#' @param haploid.genos A 2n x m matrix of haploid alleles at m loci for n
#' individuals.
#' 
#' @import dplyr
#' 
#' @export
#' 
genotypic.value <- function(genome, haploid.genos) {
  
  # Deal with inpute
  if(!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- do.call("rbind", haploid.genos)
  }
  
  # Find the number of chromosomes
  n.chr <- length(genome)
  # Find the number of loci per chromosome
  loci.per.chr <- sapply(X = genome, FUN = function(chr) chr@pos.snp %>% length())
  
  # Split the haploid.mat into chromosomes
  haploid.genos.split <- sapply(X = split(x = seq(ncol(haploid.mat)), rep(seq(n.chr), loci.per.chr)), 
                                FUN = function(chr) haploid.mat[,chr])
    
  geno.value.per.chr <- mapply(genome, haploid.genos.split, 
                               FUN = function(chr, chr.haploid.genos) {
    
    # Pull out the index of the QTL per chromosome
    qtl.index <- chr@pos.add.qtl$ID
    # Pull out the effects of the QTL
    effects <- chr@add.and.dom.eff
    add.eff <- effects$add
    dom.eff <- effects$dom
    
    # If the dominance effects are NULL, set to a vector of 0's of length n.qtl
    if (is.null(dom.eff)) { 
      dom.eff <- rep(0, length(add.eff)) 
    }
    
    # Pull out the allele states for the QTL
    qtl.alleles <- chr.haploid.genos[,qtl.index]
    
    # Split the gametes into pairs
    haploid.line.split <- split(seq(nrow(qtl.alleles)), rep(seq(nrow(qtl.alleles)/2), each = 2))
    
    # lapply to calculate genotypic value
    geno.values <- lapply(X = haploid.line.split, FUN = function(pair) {
      # Pull out the allele states for an individual
      line.gametes <- qtl.alleles[pair,]
      # Sum and subtract 1 to get the multiplier for a
      line.genos <- t(as.matrix(colSums(line.gametes) - 1))
      
      # Mulitple by the qtl allele effect to get the value at each QTL
      # Then add the domiance effect for hets
      qtl.value <- line.genos * add.eff
      qtl.value[line.genos == 0] <- dom.eff[line.genos == 0]
      
      # Sum to get the genotypic value
      sum(qtl.value)
    })
    
    # Collapse into a matrix
    geno.values <- as.matrix(do.call("rbind", geno.values))
    row.names(geno.values) <- NULL
    return(geno.values)
  })
  
  # Sum each row of the matrix
  geno.value <- rowSums(geno.value.per.chr) %>%
    as.matrix()
  return(geno.value)
  
} # Close the function