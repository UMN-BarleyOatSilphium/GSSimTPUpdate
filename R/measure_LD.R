#' Calculate linkage disequilibrium
#' 
#' @description 
#' Measures linkage disequilibrium as the correlation between QTL and markers.
#' 
#' @param genome The list of hypred genomes.
#' @param genos The n x m matrix of genotypes for n entries across m loci, 
#' including QTL.
#' @param Morgan.windwo The distance in each direction of a QTL to the markers
#' used to calculate LD. If \code{NULL}, all markers in the genome are used
#' to calculate LD with each QTL.
#' 
#' @import dplyr
#' 
#' @export
#' 
measure.LD <- function(genome,
                       genos,
                       Morgan.window = NULL
){
  
  # If a window is requested, split the genotype data by chromosome
  if (!is.null(Morgan.window)) {
    
    # Pull out the number of chromosomes
    n.chr <- length(genome)
    # Pull out the number of loci per chromsosome
    loci.per.chr <- sapply(X = genome, FUN = function(chr) length(chr@pos.snp))
    
    # Split the genotypes into chromosomes
    genos.split <- sapply(X = split(x = seq(ncol(genos)), rep(seq(n.chr), loci.per.chr)), FUN = function(i) genos[,i])
    
    # Apply a function over the genome and genotype data simulaneously
    genome.LD <- mapply(genome, genos.split, FUN = function(chr, chr.genos) {
      
      # Pull out QTL positions and index
      pos.qtl <- chr@pos.add.qtl
      # Number of loci
      n.loci <- chr@pos.snp %>%
        length()
      # Calculate the position of markers
      pos.markers <- setdiff(seq(n.loci), pos.qtl$ID) %>%
        list(ID = ., M = chr@pos.snp[.])
      
      # Find the polymorphic qtl and markers
      pos.poly.qtl <- lapply(X = pos.qtl, FUN = function(q) 
        q[apply(X = chr.genos[,pos.qtl$ID], MARGIN = 2, FUN = GSsim.TPUpdate:::is.polymorphic )] )
      pos.poly.markers <- lapply(X = pos.markers, FUN = function(m) 
        m[apply(X = chr.genos[,pos.markers$ID], MARGIN = 2, FUN = GSsim.TPUpdate:::is.polymorphic)] )
      
      # Exit if there are not polymorphic QTL or markers
      if(length(pos.poly.qtl$ID) == 0) return(NA)
      if(length(pos.poly.markers$ID) == 0) return(NA)
      
      # Iterate over the polymorphic qtl indices
      chr.LD.i <- mapply(pos.poly.qtl$ID, pos.poly.qtl$M, FUN = function(q.ID, q.M) {
        
        # Add and subtract the Morgan window, while rounding
        M.i.lower <- GSsim.TPUpdate:::round.limit(x = q.M - Morgan.window, 0, "lower")
        M.i.upper <- GSsim.TPUpdate:::round.limit(x = q.M + Morgan.window, chr@len.chr, "upper")
        # Find snps within the window
        markers.in.window <- pos.poly.markers$ID[findInterval(x = pos.poly.markers$M, vec = c(M.i.lower, M.i.upper)) == 1]
        
        # If no polymorphic markers are in the window, return NA
        if (length(markers.in.window) == 0) return(NA)
        
        # Extract geno data for those markers
        genos.markers.in.window <- chr.genos[,markers.in.window]
        
        # Correlation matrix
        cor(chr.genos[,q.ID], genos.markers.in.window)
        
      }); names(chr.LD.i) <- colnames(chr.genos)[pos.poly.qtl$ID]
      
      # Remove NAs
      chr.LD.i[!is.na(chr.LD.i)]
      
    })
    
    # Don't remove chromosome NAs
    return(genome.LD)
    
    # If the Morgan window is null, find the LD of QTL with all other genomic
    ## markers
  } else {
    
    # Pull out all QTL positions and index
    pos <- GSsim.TPUpdate:::find.pos(genome = genome)
    # QTL positions
    pos.qtl <- pos$pos.qtl
    pos.snp <- pos$pos.snp
    
    # Find the polymorphic qtl and markers
    pos.poly.qtl <- pos.qtl[apply(X = genos[,pos.qtl], MARGIN = 2, FUN = GSsim.TPUpdate:::is.polymorphic)]
    pos.poly.markers <- pos.snp[apply(X = genos[,pos.snp], MARGIN = 2, FUN = GSsim.TPUpdate:::is.polymorphic)]
    
    # Exit if nothing is polymorphic
    if(length(pos.poly.qtl) == 0) return(NA)
    if(length(pos.poly.markers) == 0) return(NA)
    
    # Allele calls of polymorphic QTL and markers
    genos.poly.qtl <- genos[,pos.poly.qtl]
    genos.poly.markers <- genos[,pos.poly.markers]
    
    # Calculate an LD matrix of QTL x markers
    LD.mat <- cor(genos.poly.qtl, genos.poly.markers)
    
    # Return the LD matrix
    return(LD.mat)
    
  } # Close the if statement
  
} # Close the function