#' Genotype loci
#' 
#' @description 
#' Creates a loci design matrix given a haploid matrix or a list of haploid 
#' matricies. This function essentially "genotypes" the individuals with
#' 100% accuracy.
#' 
#' @param genome The hypred genome.
#' @param haploid.genos A 2n x m matrix of haploid alleles at m loci for n
#' individuals.
#' @param include.QTL Logical whether to include the QTL in the design matrix.
#' Default is \code{FALSE}.
#' @param error.rate The rate of genotyping errors for a particular individual.
#' Used to draw from a binomial distribution whether to introduce error. The 
#' conversion of a heterozygous genotype to one of the homozygotes or from a 
#' homozygous genotype to a heterozygote or one of the homozygotes will occur
#' will equal probability.
#' @param missing.prop The proportion of missing data across the genotype calls.
#' Genotype calls are randomly assigned missing to obtain the desired missingness
#' proportion.
#' 
#' @import dplyr
#' @import stringr
#' 
#' @export 
#' 
genotype.loci <- function(haploid.genos, genome, include.QTL = FALSE, 
                          error.rate = 0, missing.prop = 0) {
  
  # Deal with input
  if (!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- as.matrix(do.call("rbind", haploid.genos))
  }
  
  ## Generate errors, if the rate is greater than 0
  if (error.rate > 0) {
    n.loci <- ncol(haploid.mat)
    n.entries <- nrow(haploid.mat)
    
    # Designate calls for errors
    error.assign <- rbinom(n = n.calls, size = 1, prob = error.rate) %>%
      matrix(nrow = nrow(haploid.mat), ncol = ncol(haploid.mat), byrow = T)
    
    haploid.mat1 <- haploid.mat - error.assign %>%
      abs()
  
  } else { # Otherwise just rename the matrix
    haploid.mat1 <- haploid.mat
  }
  
  ## Convert to genotypes
  # Pull out number of chromosomes
  n.chr <- length(genome)
  # Pull out the number of loci per chromosome
  loci.per.chr <- sapply(X = genome, FUN = function(chr) chr@pos.snp %>% length())
  # Pull out the position of the QTL on each chromosome
  qtl.index.per.chr <- sapply(X = genome, FUN = function(chr) chr@pos.add.qtl$ID )
  
  # Pull out line names
  line.names <- row.names(haploid.mat1) %>%
    # Remove gamete indicator
    str_replace(pattern = "\\.[0-9]$", replacement = "") %>%
    unique()
  
  # Split the haploid.mat into chromosomes
  haploid.genos.split <- sapply(X = split(x = seq(ncol(haploid.mat1)), rep(seq(n.chr), loci.per.chr)), 
                                FUN = function(chr) haploid.mat1[,chr])
  
  if (!include.QTL) {
    
    # Remove QTL from the haploid matrix
    haploid.genos.no.qtl <- mapply(haploid.genos.split, qtl.index.per.chr, 
                                   FUN = function(chr.haploid.genos, chr.qtl.ind)
                                     chr.haploid.genos[,-chr.qtl.ind] )

    # Reform the haploid matrix
    haploid.mat2 <- do.call("cbind", haploid.genos.no.qtl)
    
  } else { # If the genotypes of the QTL are desired
    
    haploid.mat2 <- haploid.mat1
    
  } # Close the QTL inclusion if statement
  
  # Split the haploid matrix into pairs of rows
  haploid.line.split <- split(seq(nrow(haploid.mat2)), rep(seq(nrow(haploid.mat2)/2), each = 2))
  
  # Apply a function to generate -1, 0, 1 genotype values for each pair of rows
  geno.mat <- lapply(X = haploid.line.split, FUN = function(pair) {
    # Subset the gamete matrix
    line.haploids <- haploid.mat2[pair,]
    # Find the number of 1 alleles
    line.alleles <- colSums(line.haploids)
    # Subtract one to get the coded genotypes
    line.alleles - 1
  }) %>%
    do.call("rbind", .)
  
  ## Add missing data, if the proportion is greater than 0
  if (missing.prop > 0) {
    
    n.calls = nrow(geno.mat) * ncol(geno.mat)
    
    # Designate calls for missing
    missing.assign <- rbinom(n = n.calls, size = 1, prob = missing.prop) %>%
      as.logical() %>%
      matrix(nrow = nrow(haploid.mat), ncol = ncol(haploid.mat), byrow = T)
    
    # Set those calls as missing
    geno.mat1 <- geno.mat
    geno.mat1[missing.assign] <- NA
    # Refit to matrix
    geno.mat2 <- geno.mat1 %>%
      matrix(nrow = nrow(geno.mat), ncol = ncol(geno.mat))
    
  } else {
    geno.mat2 <- geno.mat
  }
  
  # Add line names to the rows
  row.names(geno.mat2) <- line.names
  colnames(geno.mat2) <- colnames(geno.mat)
    
  # Return the matrix
  return(geno.mat2)
  
} # Close the function