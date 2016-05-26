# Simulation functions source code
# Jeff Neyhart
# March 15, 2016

###### Define functions #####

# Define a function to convert a hmp format of genotypes to a genotype matrix
hmp2mat <- function(hmp) {
  mat <- t(hmp[,-c(1:4)])
  dimnames(mat) <- list(colnames(hmp)[-c(1:4)], hmp[,1])
  return(as.matrix(mat))
} # Close the function

  

# Define a function to calculate pairwise per-SNP diversity from a genotype matrix
# This is taken from Tom Kono
# By default, the script will remove missing values and heterozygotes
SNP.pairwise.div <- function(geno.mat, 
                             missing = c(NA, 0)) {
  
  # Apply a function over the genotype matrix
  per.SNP.div <- apply(X = geno.mat, MARGIN = 2, function(x) {
    #   First, remove missing data
    #   is.element(a, b) will return a vector of boolean values of length
    #   length(a), indicating whether or not a[i] is in b.
    #   This removes heterozygous calls and missing calls from our genotypes.
    calls <- x[!is.element(x, missing)]
    #   Count up the genotypic classes
    counts <- table(calls)
    
    #   If the length of the counts is 1, then the SNP is monomorphic, and
    #   pairwise diversity is 0.
    if(length(counts) == 1) {
      pairwise_div <- 0
    } else {
      #   Then get the minimum and max for major and minor alleles
      major <- max(counts)
      minor <- min(counts)
      #   The way to calculate pairwise *similarity* for a single SNP is
      #   [(MinorCount choose 2) + (MajorCount choose 2)] / (N.Ind choose 2)
      #   Essentially, the number of similarities due to major/major
      #   comparisons and minor/minor comparisons, divided by the total number
      #   of comparisons.
      minorsim <- choose(minor, 2)
      majorsim <- choose(major, 2)
      nindsim <- choose(length(calls), 2)
      pairwise_sim <- (minorsim + majorsim) / nindsim
      #   Pairwise div. is 1-similarity
      pairwise_div <- 1 - pairwise_sim
    }
    # Return the diversity measure
    return(pairwise_div)
  })
  
  # Return the diversity across all SNPs
  return(per.SNP.div)
} # Close the function

# Define a function to create a family
## This function requires the genome, a 2 x m matrix of the parent 1 genome, a 2 x m matrix of the parent 2 genome,
## the size of the family, and mutation rates
make.family <- function(genome, # hypred genome 
                        parent1.genome, # 2 x m matric of parent 1 genome
                        parent2.genome, # 2 x m matrix of parent 2 genome
                        N, # Number of individuals in the family
                        generations, # Number of generations to advance the family beyond the F_1. For an inbred population, the result will be a F_n family, where n = generations + 1
                        pop.type = "random", # The type of population. By default it is random
                        cycle.number,
                        family.number,
                        mutate = TRUE,
                        mutation.rate.snp, 
                        mutation.rate.qtl) {
  
  
  
  # Determine the number of gametes to produce
  n.gamete <- N * 2
  
  F_1 <- do.call("rbind", lapply(X = 1:N, function(ind) {
    # Gamete 1
    gamete1 <- recombine(genome = genome,
                         haploid.genomeA = parent1.genome[1,],
                         haploid.genomeB = parent1.genome[2,],
                         mutate = mutate,
                         mutation.rate.snp = mutation.rate.snp,
                         mutation.rate.qtl = mutation.rate.qtl)
    
    # Gamete 2
    gamete2 <- recombine(genome = genome,
                         haploid.genomeA = parent2.genome[1,],
                         haploid.genomeB = parent2.genome[2,],
                         mutate = mutate,
                         mutation.rate.snp = mutation.rate.snp,
                         mutation.rate.qtl = mutation.rate.qtl)
    
    return(rbind(gamete1, gamete2)) }) )
  
  # Iterate changes to the population over a number of generations
  # Rename the population
  pop.i <- F_1
  
  # If no advancement is requested, don't advance
  if (generations > 0) {
    
    # Iterate over generations
    for (g in 1:generations) {
      
      # This applies the function over individuals
      # First we generate a set of recombinant gametes
      recomb.pop <- lapply(X = seq(1, n.gamete, 2), FUN = function(i) {
        gamete.index1 <- i
        gamete.index2 <- i + 1
        
        # Simulate the gametes from the i-th individual
        t(replicate(2, recombine(genome = genome,
                                 haploid.genomeA = pop.i[gamete.index1,],
                                 haploid.genomeB = pop.i[gamete.index2,],
                                 mutate = mutate,
                                 mutation.rate.snp = mutation.rate.snp,
                                 mutation.rate.qtl = mutation.rate.qtl))) })
        
      # Collapse the list into the correct matrix
      pop.i <- do.call("rbind", recomb.pop)
      
      # If randomly mating, permutate
      if (pop.type == "random") {
        pop.i <- pop.i[sample(1:n.gamete),]
      }
      
    } # Close the generation for loop
    
  } # Close the if statement
  
  # Name the individuals
  row.names(pop.i) <- paste("C", cycle.number, "_", (1+generations), formatC(family.number, width = 3, format = "d", flag = "0"), "-", formatC(rep(c(1:N), each = 2), width = 4, format = "d", flag = "0"), ".", rep(c(1:2), 2), sep = "")
  
  # Return the population matrix
  return(pop.i)
  
} # Close the function

# Define a function to advance or change a family
## This function requires a genome, the 2N x m population matrix, the type of family ("random" or "inbred"),
## the number of generations, and mutation rates
advance.family <- function(genome, 
                           pop.mat, 
                           pop.type = "random",
                           cycle.number,
                           starting.generation, 
                           generations, 
                           mutation.rate.snp, 
                           mutation.rate.qtl) {
  
  # Deal with inputs
  pop.mat <- as.matrix(pop.mat)
  # Gather marker names
  markers <- colnames(pop.mat)
  
  # Rename the population
  pop.i <- pop.mat
  
  # Split line information on the inbreeding generation
  line.names <- row.names(pop.i)
  line.name.length <- unique(sapply(X = line.names, FUN = nchar))
  
  prefix <- paste("C", cycle.number, sep = "")
  suffix <- substring(text = line.names, first = (line.name.length - 9))
  
  # Find the number of gametes
  n.gamete <- nrow(pop.i)
  
  # If no advancement is requested, don't advance
  if (generations > 0) {
    # Iterate over generations
    for (g in 1:generations) {
      
      # This applies the function over individuals
      # First we generate a set of recombinant gametes
      recomb.pop <- lapply(X = seq(1, n.gamete, 2), FUN = function(i) {
        gamete.index1 <- i
        gamete.index2 <- i + 1
        
        # Simulate the gametes from the i-th individual
        gamete1 <- 
          hypredRecombine(genome,
                          genomeA = pop.i[gamete.index1,],
                          genomeB = pop.i[gamete.index2,],
                          mutate = T,
                          mutation.rate.snp = mutation.rate.snp,
                          mutation.rate.qtl = mutation.rate.qtl,
                          block = FALSE)
        
        gamete2 <-
          hypredRecombine(genome,
                          genomeA = pop.i[gamete.index1,],
                          genomeB = pop.i[gamete.index2,],
                          mutate = T,
                          mutation.rate.snp = mutation.rate.snp,
                          mutation.rate.qtl = mutation.rate.qtl,
                          block = FALSE)
        
        return( list(gamete1, gamete2) )
      })
      
      # Collapse the list into the correct matrix
      pop.i <- do.call("rbind", lapply(X = recomb.pop, FUN = function(x) do.call("rbind", x)))
      
      # If randomly mating, permutate
      if (pop.type == "random") {
        pop.i <- pop.i[sample(1:n.gamete),]
      }
      
    } # Close the generation for loop
  } # Close the if statement
  
  # Assign line names
  row.names(pop.i) <- paste(prefix, "_", (starting.generation + generations), suffix, sep = "")
  # Assign marker names
  colnames(pop.i) <- markers
  
  # Return the population matrix
  return(pop.i)
} # Close the function

# Define a function to convert observed genotypes to gamete codes
geno.reverse.code <- function(geno.mat) {
  
  # Deal with inpute
  geno.mat <- as.matrix(geno.mat)
  
  # Pull out line names
  line.names <- row.names(geno.mat)
  
  # Find the number of entried and markers
  n <- nrow(geno.mat)
  m <- ncol(geno.mat)
  
  gamete.mat <- apply(X = geno.mat, MARGIN = 1, FUN = function(geno) {
    # Add 1 to get the number of primary alleles at each locus
    geno <- round(geno) + 1
    # Apply a function over the vector of genotypes to create gametes
    gamete_i <- sapply(X = geno, FUN = function(g_i) {
      if(g_i == 0) return(c(0,0))
      if(g_i == 1) return(c(0,1))
      if(g_i == 2) return(c(1,1))
    })
    # Return the gametes
    return(list(gamete_i))
  })
  
  # Collapse the list into a matrix
  gamete.mat <- do.call("rbind", lapply(X = gamete.mat, FUN = function(i) return(i[[1]])))
  
  # Add line names back in
  line.names <- paste(rep(line.names, each = 2), ".", c(1,2), sep = "")
  row.names(gamete.mat) <- line.names
  
  # Return matrix
  return(gamete.mat)
} # Close the function

# Define a function to change heterozygous markers to homozygous based on the flanking markers
# The function will convert hets (0) to either 1 or -1, and will weight the probability of each genotype
## on the distance from the current marker each flanking marker
# The function requires a n x m genotype matrix and a m x 3 data.frame of marker name, chromosome, and cM position
het.impute <- function(geno.mat, 
                       marker.map) {
  
  # Apply the function across rows (genotypes)
  # This function will return a genotype matrix without hets
  geno.sin.hets <- t(apply(X = geno.mat, MARGIN = 1, FUN = function(geno) {
    # Find the positions of loci with a het call
    het.index <- which(geno == 0)
    
    # While loop to continue as long as there are hets
    while(length(het.index) > 0) {
      # Apply a new function over the indices
      # This function will return a vector of the het markers and their "imputed" genotypes
      het.impute <- sapply(X = het.index, FUN = function(l) {
        # Extract the marker name and position
        het.marker <- colnames(geno.mat)[l]
        het.marker.pos <- marker.map[marker.map[,1] %in% het.marker,c(2,3)]
        
        # Extract the marker names of the flanking markers
        if (l - 1 > 0) {
          left.marker <- colnames(geno.mat)[l - 1]
          left.marker.pos <- marker.map[marker.map[,1] %in% left.marker,c(2,3)]
          # Find the cM distance between the left flanking markers
          # Only do this if the chromosomes are the same
          if (left.marker.pos[1] == het.marker.pos[1]) {
            left.marker.dist <- as.numeric(abs(left.marker.pos[2] - het.marker.pos[2]))
            # Assign a probabilty to the marker based on its distance
            left.marker.prob <- 1 - (left.marker.dist / 100)
            # Pull out the genotype for the left marker
            left.marker.geno <- as.numeric(geno[l - 1])
            if (left.marker.geno == 0) {left.marker.prob <- 0}
          } else {
            left.marker.prob <- 0
            left.marker.geno <- NA
          }
          # If the "left marker" does not exist, set it's probability to 0
        } else {
          left.marker.prob <- 0
          left.marker.geno <- NA
        }
        
        # Repeat the if statement for the right marker
        if (l + 1 < length(geno)) {
          right.marker <- colnames(geno.mat)[l + 1]
          right.marker.pos <- marker.map[marker.map[,1] %in% right.marker,c(2,3)]
          # Find the cM distance between the right flanking markers
          # Only do this if the chromosomes are the same
          if (right.marker.pos[1] == het.marker.pos[1]) {
            right.marker.dist <- as.numeric(abs(right.marker.pos[2] - het.marker.pos[2]))
            # Assign a probabilty to the marker based on its distance
            right.marker.prob <- 1 - (right.marker.dist / 100)
            # Pull out the genotype for the right marker
            right.marker.geno <- as.numeric(geno[l + 1])
            # If the flanking marker is also a het, give a probability of 0
            if (right.marker.geno == 0) {right.marker.prob <- 0}
          } else {
            right.marker.prob <- 0
            right.marker.geno <- NA
          }
          # If the "right marker" does not exist, set it's probability to 0
        } else {
          right.marker.prob <- 0
          right.marker.geno <- NA
        }
        
        # If the flanking probabilities are both 0, set the genotype as a het
        if(all(right.marker.prob == 0, left.marker.prob == 0)) {
          het.marker.geno <- 0
        } else {
          # Assign the het.marker genotype by randomly sampling the flanking marker genotypes
          ## and using a probability of each genotype that is proportional to their distance from 
          ## the het marker
          het.marker.geno <- sample(x = c(left.marker.geno, right.marker.geno), size = 1, prob = c(left.marker.prob, right.marker.prob) )
        }
        
        # Return the genotype
        return(het.marker.geno)
      })
      
      # Reassign the het markers according to the imputed genotypes
      geno[het.index] <- het.impute
      # Check again for hets
      het.index <- which(geno == 0)
    } # Close the while loop
    
    # Return the "imputed" genos
    return(geno)
  })) # Close the apply function
  
  # Return the new matrix
  return(geno.sin.hets)
} # Close the function

# Define a function to create a crossing block based on a vector of line names for each parents (this can be the same vector) 
# and a number of crosses, while choosing whether or not to allow reciprocal crosses
make.crossing.block <- function(parent1.lines, # Character vector of lines for the first parent
                                parent2.lines, # Character vector of lines for the second parent
                                n.crosses, # Number of crosses to make
                                method = "random", # Method to assign parents. Can be "random" for random pairs of the parent1 and parent2, "chain" for sequential crosses (Note: can only be used if the parent1.line and parent2.lines vectors are identical), or "all.pairwise" for all possible pairwise crosses.
                                use.parents.once = FALSE) {
  
  # Find the intersection of parent1.lines and parent2.lines
  par.intersect <- intersect(parent1.lines, parent2.lines)
  n.intersect <- length(par.intersect)
  # Find the number of parent1 and parent2 lines
  n.par1 <- length(parent1.lines)
  n.par2 <- length(parent2.lines)
  
  if ( all(n.intersect == n.par1, n.intersect == n.par2) ) {
    n.possible.crosses <- (n.par1 * (n.par2 - 1)) / 2
    same.lines = T
  } else {
    n.possible.crosses <- n.par1 * n.par2
    same.lines = F
  }
  
  # Create all pairwise crosses
  sample.crosses <- expand.grid(parent1.lines, parent2.lines)[,c(2,1)]
  # Remove selfs
  sample.crosses <- subset(x = sample.crosses, !apply(X = sample.crosses, MARGIN = 1, FUN = function(cross) any(duplicated(cross))) )

  # If statements for methods
  if (method == "all.pairwise") {
    return(sample.crosses)
  }
  if (method == "random") {
    # First see if the number of requested crosses is more than possible
    if (n.crosses > n.possible.crosses) stop("n.crosses is more than possible.")
    
    # If parents should only be used once, sample without replacement all lines and put them into a matrix
    if (use.parents.once) {
      # Sample into a matrix
      if (same.lines) {
        random.crosses <- as.data.frame(matrix(sample(parent1.lines), ncol = 2, byrow = T))
        crosses.ind <- sort(sample(1:nrow(random.crosses), n.crosses))
        random.crosses <- random.crosses[crosses.ind,]
      } else {
        # Are there any overlapping lines in the two parent vectors?
        if (n.intersect == 0) {
          # If so just sample each vector into a matrix
          random.crosses <- as.data.frame(cbind( sample(parent1.lines), sample(parent2.lines) ))
          crosses.ind <- sort(sample(1:nrow(random.crosses), n.crosses))
          random.crosses <- random.crosses[crosses.ind,]
        } else {
          stop("Function not available.")
        }
      }
      
    } else {
      
      # Remove reciprocal
      sample.crosses <- subset(x = sample.crosses, subset = !duplicated(t(apply(sample.crosses, MARGIN = 1, FUN = sort))))
      crosses.ind <- sort(sample(1:nrow(sample.crosses), n.crosses))
      random.crosses <- sample.crosses[crosses.ind,]
    }
    
    colnames(random.crosses) <- c("Parent1", "Parent2")
    row.names(random.crosses) <- NULL
    return(random.crosses)
  }
  
} # Close the function

# Define a function to take gamete information and a crossing block and produce a list families
# The function also requires the genome, the size of the families and the mutation rates
make.population <- function(genome, 
                            parental.haploids, 
                            crossing.block, 
                            N, 
                            cycle.number, 
                            generations, 
                            pop.type = "random", 
                            mutation.rate.snp, 
                            mutation.rate.qtl) {
  
  # Make sure the crossing block is a matrix
  crossing.block <- as.matrix(crossing.block)
  # Pull out the loci names
  loci.names <- colnames(parental.haploids)
  
  # Apply a function over the crossing block
  # This function will return a list of the cross families
  fam.list <- lapply(X = 1:nrow(crossing.block), FUN = function(i) {
    
    # Identify the cross
    cross <- crossing.block[i,]
    
    # Use grep to pull out the gametes
    parent1.haploid <- parental.haploids[grep(pattern = cross[1], x = row.names(parental.haploids)),]
    parent2.haploid <- parental.haploids[grep(pattern = cross[2], x = row.names(parental.haploids)),]
    
    # Make a population from each cross
    pop <- make.family(genome = genome, 
                       parent1.genome = parent1.haploid, 
                       parent2.genome = parent2.haploid, 
                       N = N, 
                       generations = generations, 
                       pop.type = pop.type,
                       family.number = i,
                       cycle.number = cycle.number,
                       mutation.rate.snp = mutation.rate.snp, 
                       mutation.rate.qtl = mutation.rate.qtl)
    
    # Add marker names back in
    colnames(pop) <- loci.names
    
    # Return the population genotypes in a list
    return(pop)
  })
  
  # Rename the list entries
  names(fam.list) <- apply(X = crossing.block, MARGIN = 1, FUN = function(cross) return(paste(cross[1],cross[2], sep = ".")))
  
  # Return the list
  return(fam.list)
} # Close the function

# Define a function to take a list of populations from a cross and advance each population either through
## random mating or through inbreeding
advance.population <- function(genome, 
                               list.of.families, 
                               pop.type = "random", 
                               generations,
                               cycle.number,
                               starting.generation,
                               mutation.rate.snp, 
                               mutation.rate.qtl) {
  
  # Apply a function over the list
  advanced.family <- lapply(X = list.of.families, FUN = function(family) {
    advanced.family <- advance.family(genome = genome, 
                                      pop.mat = family, 
                                      pop.type = pop.type, 
                                      generations = generations, 
                                      cycle.number = cycle.number, 
                                      starting.generation = starting.generation,
                                      mutation.rate.snp = mutation.rate.snp,
                                      mutation.rate.qtl = mutation.rate.qtl)
    return(advanced.family)
  })
  
  return(advanced.family)
} # Close the function

# Define a function to return the genotypic value of a set of lines based on their QTL alleles and the effect of those alleles
genotypic.value <- function(genome, 
                            haploid.genos) {
  
  # Deal with inpute
  if(!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- do.call("rbind", haploid.genos)
  }
  
  # Find the number of chromosomes
  n.chr <- length(genome)
  # Find the number of loci per chromosome
  loci.per.chr <- sapply(X = genome, FUN = function(chromosome) length(slot(chromosome, "pos.snp")))
  
  # Split the haploid.mat into chromosomes
  haploid.chr.split <- split(x = 1:ncol(haploid.mat), rep(1:n.chr, loci.per.chr))
  haploid.list <- sapply(X = haploid.chr.split, FUN = function(chr) haploid.mat[,chr] )
  
  # Apply a function over each chromosome in the genome
  geno.value.per.chr <- sapply(X = 1:n.chr, FUN = function(i) {
    
    # Pull out the chromosome
    chromosome <- genome[[i]]
  
    # Pull out the index of the QTL per chromosome
    qtl.index <- slot(chromosome, "pos.add.qtl")$ID
    # Pull out the effects of the QTL
    effects <- slot(chromosome, "add.and.dom.eff")
    add.eff <- effects$add
    dom.eff <- effects$dom
  
    # If the dominance effects are NULL, set to a vector of 0's of length n.qtl
    if (is.null(dom.eff)) { 
      dom.eff <- rep(0, length(add.eff)) 
    }
    
    # Pull out the allele states for the QTL
    qtl.alleles <- haploid.list[[i]][,qtl.index]
  
    # Split the gametes into pairs
    haploid.line.split <- split(1:nrow(qtl.alleles), rep(1:(nrow(qtl.alleles)/2), each = 2))
    
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
  geno.value <- as.matrix(rowSums(geno.value.per.chr))
  return(geno.value)
  
} # Close the function
  
  
  

# Define a function to generate value measures on the population of interest. The values can be the genotypic values (sum
# of the a values at each locus) or the phenoypic values (genetic values + error)
# This function will require the genome, the matrix, heritability, and distribution of phenotypes
evaluate.population <- function(genome, 
                                haploid.genos, 
                                h2,
                                just.geno.values = F, # Should the function only output the genotypic values?
                                V_e.scale = 8, # By how much should the environmental variance be larger than the genotypic variance?
                                n.env = 3,
                                n.rep = 2
) {
  
  # Deal with input
  if(!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- do.call("rbind", haploid.genos)
  }
  
  # Pull out the line names
  line.names <- row.names(haploid.mat)
  # Condense the names
  line.names <- unique(sub(pattern = "\\.[0-9]$", replacement = "", x = line.names))
  
  # First generate genotype values based on the genome and genotypes
  g <- genotypic.value(genome = genome, haploid.genos = haploid.mat)
  row.names(g) <- line.names
  # Find the length (i.e. number of genotypes)
  n.geno = length(g)
  
  # Calculate the true genetic variance
  V_g <- var(g)
  
  # Determine the residual variance from the heritability
  # V_residual is V_ge and V_epsilon confounded
  V_residual = n.rep * n.env * ((V_g / h2) - V_g)
  
  # If just genotypic values are requested, stop and return a list
  if (just.geno.values) {
    output.list <- list(geno.values = g,
                        var.components = list(V_g = V_g, V_residual = V_residual))
    return(output.list)
  }
  
  # Calculate the environmental variance as a multiplicative scale of the genetic variance
  V_e <- V_g * V_e.scale
  # Sample the environmental effects
  e = rnorm(n.env, 0, sqrt(V_e))
  
  # Sample the residuals
  # Remember this includes error and gxe effect
  residual <- matrix(rnorm(n.geno * n.env * n.rep, 0, sqrt(V_residual)), n.geno, n.env * n.rep)
  
  # Calculate phenotypic values
  y <- matrix(g, n.geno, (n.env * n.rep)) + matrix(e, n.geno, (n.env * n.rep), byrow = T) + residual
  
  # Label the matrix
  row.names(y) <- line.names
  colnames(y) <- paste( paste("env", 1:n.env, sep = ""), rep(paste("rep", 1:n.rep, sep = ""), each = n.env), sep = "." )
  
  # Find the average across environments as the value to use
  mu_y <- as.matrix(rowMeans(y))
  row.names(mu_y) <- line.names
  
  # Find the average across genotypes
  mu.p <- mean(mu_y) # Phenotypic value
  mu.g <- mean(g) # Genotypic values
  
  # Output list
  output.list <- list(full.matrix = y,
                      geno.values = g,
                      mean.pheno.values = mu_y,
                      mu.p = mu.p,
                      mu.g = mu.g,
                      var.components = list(V_g = V_g, V_e = V_e, V_residual = V_residual))
  
  return(output.list)
  
} # Close the function


# Define a function to subset a list of values
subset.values <- function(values.list,
                          lines.to.subset) {
  
  # Subset pheno and geno values
  subset.pheno.values <- as.matrix(values.list$mean.pheno.values[lines.to.subset,])
  subset.geno.values <- as.matrix(values.list$geno.values[lines.to.subset,])
  
  # Calculate mean pheno and geno values
  subset.mu.p <- mean(subset.pheno.values)
  subset.mu.g <- mean(subset.geno.values)
  
  # Calculate genetic variance of the subset
  subset.V_g <- var(subset.geno.values)
  
  # Package and return a list
  subset.list <- list(geno.values = subset.geno.values,
                      mean.pheno.values = subset.pheno.values,
                      mu.p = subset.mu.p,
                      mu.g = subset.mu.g,
                      var.components = list(true = list(V_g = subset.V_g)))
  return(subset.list)
  
} # Close the function
  


# Define a function to apply selection
# This function will require a vector of phenotypes and a selection intensity (proportion of samples to select)
# The selection can be "best", "worst", or "random"
select.population <- function(pheno.mat, 
                              sel.intensity,
                              selection = "best",
                              exclusion = NULL # A vector of line names to exclude when making selections
                              ) {
  
  # Exclude based on the provided exclusion vector
  pheno.mat <- as.matrix(pheno.mat[!row.names(pheno.mat) %in% exclusion,])
  
  # Number in the pheno.vec
  N.tot <- nrow(pheno.mat)
  
  # If the sel.intensity is between 0 and 1, find the total number to keep
  if (all(sel.intensity > 0, sel.intensity < 1)) {
    # Number to select
    N.sel <- sel.intensity * N.tot
  } else {
    N.sel <- sel.intensity
  }
  
  # If the selection does not have a "+" (i.e. random, best, worst), continue appropriately
  if (!grepl(pattern = "\\+", x = selection)) {
    # If using the "best" method, select the top N.sel
    if (selection == "best") {
      # Sort the phenotypes
      pheno.mat.ordered <- as.matrix(pheno.mat[order(pheno.mat, decreasing = T),])
      pheno.sel <- pheno.mat.ordered[1:N.sel,]
    }
    # If using the "worst" method, select the bottom N.sel
    if (selection == "worst") {
      # Sort the values
      pheno.mat.ordered <- as.matrix(pheno.mat[order(pheno.mat, decreasing = F),])
      pheno.sel <- pheno.mat.ordered[1:N.sel,]
    }
    # If using "random", randomly select the lines
    if (selection == "random") {
      pheno.sel <- pheno.mat[sample(1:N.tot, size = N.sel),]
    }
  
  } else { # Otherwise return all lines
    pheno.sel <- pheno.mat[sample(1:N.tot, size = N.tot),]
  }
  
  # Sort the selections based on name
  pheno.sel <- pheno.sel[order(names(pheno.sel))]
  # Retrieve the names
  line.names <- names(pheno.sel)
  
  # Create output list
  output.list <- list(value.sel = as.matrix(pheno.sel), lines.sel = line.names)
  
  # Return the selections
  return(output.list)
} # Close the function

# Define a function to subset a gamete matrix or list of gamete matricies based on a character vector of line names
subset.gametes <- function(gametes,
                           line.names) {
  
  # Deal with inputs
  if (!is.list(gametes)) {
    gamete.mat <- as.matrix(gametes)
  } else {
    gamete.mat <- do.call("rbind", gametes)
  }
  
  # Shorten the query lines name to just the cross name
  line.names <- as.character(line.names)
  shortened.line.names <- substring(text = line.names, (nchar(line.names) - 7), nchar(line.names))
  
  # Extract names from the gamete.mat matrix
  gamete.lines <- row.names(gamete.mat)
  # Remove the gamete indicator at the end of the name
  gamete.lines <- sub(pattern = "\\.[0-9]$", replacement = "", x = gamete.lines)
  # Shorten the gamete line names to the cross name
  shortened.gamete.names <- substring(text = gamete.lines, (nchar(gamete.lines) - 7), nchar(gamete.lines))
  
  # Find an index of line.names within the gamete.lines
  line.index <- shortened.gamete.names %in% shortened.line.names
  
  # Subset the gamete matrix using the index
  gamete.subset <- subset.matrix(x = gamete.mat, subset = line.index)
  # Return the subset
  return(gamete.subset)
} # Close the function

# Define a function to create a loci design matrix given a haploid matrix or a 
## list of haploid matricies. This function essentially "genotypes" the family 
##or population for loci
genotype.loci <- function(haploid.genos,
                          genome,
                          include.QTL = FALSE) {
  
  # Deal with input
  if (!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- as.matrix(do.call("rbind", haploid.genos))
  }
  
  # Pull out number of chromosomes
  n.chr <- length(genome)
  # Pull out the number of loci per chromosome
  loci.per.chr <- sapply(X = genome, FUN = function(chromosome) length(slot(chromosome, "pos.snp")))
  # Pull out the position of the QTL on each chromosome
  qtl.index.per.chr <- sapply(X = genome, FUN = function(chromosome) slot(chromosome, "pos.add.qtl")$ID)
  
  # Pull out line names
  line.names <- unique(sub(pattern = "\\.[0-9]$", replacement = "", x = row.names(haploid.mat)))
  
  # Split the haploid.mat into chromosomes
  haploid.chr.split <- split(x = 1:ncol(haploid.mat), rep(1:n.chr, loci.per.chr))
  haploid.list <- sapply(X = haploid.chr.split, FUN = function(chr) haploid.mat[,chr] )
  
  if (!include.QTL) {
    
    # Remove QTL from the haploid matrix
    haploid.list.no.qtl <- sapply(X = 1:n.chr, FUN = function(i) haploid.list[[i]][,-qtl.index.per.chr[[i]]] )
    # Reform the haploid matrix
    haploid.mat.no.qtl <- do.call("cbind", haploid.list.no.qtl)
    
    # Split the haploid matrix into pairs of rows
    haploid.line.split <- split(1:nrow(haploid.mat.no.qtl), rep(1:(nrow(haploid.mat.no.qtl)/2), each = 2))
    
    # Apply a function to generate -1, 0, 1 genotype values for each pair of rows
    geno.list <- lapply(X = haploid.line.split, FUN = function(pair) {
      # Subset the gamete matrix
      line.haploids <- haploid.mat.no.qtl[pair,]
      # Find the number of 1 alleles
      line.alleles <- colSums(line.haploids)
      # Subtract one to get the coded genotypes
      line.alleles - 1
    })
    
    # Collapse
    geno.mat <- do.call("rbind", geno.list)
    
    # Add line names to the rows
    row.names(geno.mat) <- line.names
    
  } else { # If the genotypes of the QTL are desired
    
    # Split the haploid matrix into pairs of rows
    haploid.line.split <- split(1:nrow(haploid.mat), rep(1:(nrow(haploid.mat)/2), each = 2))
    
    # Apply a function to generate -1, 0, 1 genotype values for each pair of rows
    geno.list <- lapply(X = haploid.line.split, FUN = function(pair) {
      # Subset the gamete matrix
      line.haploids <- haploid.mat[pair,]
      # Find the number of 1 alleles
      line.alleles <- colSums(line.haploids)
      # Subtract one to get the coded genotypes
      line.alleles - 1
    })
    
    # Collapse
    geno.mat <- do.call("rbind", geno.list)
    
    # Add line names to the rows
    row.names(geno.mat) <- line.names
    
  }
  
  # Return the matrix
  return(geno.mat)
} # Close the function


# Define a function to conduct genomic predictions
# The function requires training population phenotypes, a design matrix of training marker genotypes, and design
## matrix of candidate marker genotypes.
# The output of the functtion is a matrix of GEBVs and the solution list to the mixed model
make.predictions <- function(pheno.train,
                             geno.train,
                             geno.pred,
                             model = "RRBLUP") {
  
  # Deal with input
  pheno.train <- as.matrix(pheno.train)
  geno.train <- as.matrix(geno.train)
  geno.pred <- as.matrix(geno.pred)
  
  # Find number of lines and markers
  n.TP <- nrow(geno.train)
  n.CP <- nrow(geno.pred)
  m <- ncol(geno.train)
  
  # Create model parameters
  y <- as.numeric(pheno.train)
  X <- matrix(1, nrow = n.TP)
  # Create model parameters based on requested model
  if (model == "RRBLUP") {
    Z <- geno.train
    K <- diag(m)

    # Solve the mixed model
    system.time(solve.out <- mixed.solve(y = y, X = X, Z = Z))
    marker.effects <- solve.out$u
    # Calculate GEBVs
    GEBV <- geno.pred %*% marker.effects
  } 
  
  if (model == "GBLUP") {
    Z <- diag(n.TP + n.CP)
    Z <- as.matrix(Z[1:n.TP,])
    K <- A.mat(X = rbind(geno.train, geno.pred), 0, 1)
    
    # Solve the mixed model
    system.time(solve.out <- mixed.solve(y = y, X = X, Z = Z, K = K))
    GEBV <- solve.out$u[-c(1:n.TP)]
  }
  
  # Convert GEBVs to matrix
  GEBV <- as.matrix(GEBV)
  row.names(GEBV) <- row.names(geno.pred)
  
  # Create an output list
  output.list <- list(GEBV = GEBV, 
                      matricies = list(y = y, X = X, Z = Z, K = K), 
                      solve.out = solve.out, 
                      parameters = list(n.TP = n.TP, n.CP = n.CP, n.markers = m))
  # Return the data
  return(output.list)
} # Close the function


# Define a function to validate predictions
# The function requires training and validation population phenotypes, and design matrices of training and
## validation marker genotypes
validate.predictions <- function(predicted.GEBVs,
                                 observed.values,
                                 boot.reps = NULL) {
  
  # Deal with input
  predicted.GEBVs <- as.matrix(predicted.GEBVs)
  observed.values <- as.matrix(observed.values)
  
  # If no boot.reps, do not perform bootstrapping
  if (is.null(boot.reps)) {
    
    pred.r <- cor(predicted.GEBVs, observed.values)
    pred.r.sd <- NULL
    
  } else { # Otherwise perform bootstrapping
  
    # Perform bootstrapping to measure the correlation between GEBVs and phenotypes
    # First combine the data
    data.to.bootstrap <- cbind(predicted.GEBVs, observed.values)
    # Define a function for performing a correlation of the sample given by the ith replication
    boot.cor <- function(input.data, i) {
      rep.data <- input.data[i,]
      rep.cor <- cor(rep.data[,1], rep.data[,2])
      return(rep.cor)
    }
    # Perform bootstrapping
    boot.results <- boot(data = data.to.bootstrap, statistic = boot.cor, R = boot.reps)
    
    # Parse the results
    pred.r <- boot.results$t0
    pred.r.sd <- sd(boot.results$t)
  }
  
  # Create an output list
  output.list <- list(pred.r = pred.r, pred.r.sd = pred.r.sd, boot.reps = boot.reps)
  # Return the data
  return(output.list)
} # Close the function


# Define a function to calculate allele frequencies
# The function will calculate the frequency of the 1 allele across all genotypes in a genotype matrix
# It will also plot the site-frequency spectrum, if desired
calculate.allele.freq <- function(geno.mat = NULL,
                                  plot.sfs = FALSE) {
  
  # Deal with input
  geno.mat <- as.matrix(geno.mat)
  
  # Add 1 to get the count of the first allele at each locus per entry
  geno.mat <- geno.mat + 1
  
  # Apply a function over the columns of the matrix (SNPs)
  freq.alleles <- apply(X = geno.mat, MARGIN = 2, FUN = function(snp) {
    # Find the total number of "1" alleles
    sum.1 <- sum(snp)
    # Find the total number of alleles
    n.alleles <- 2 * length(snp)
    # Calculate frequency
    freq.1 <- sum.1 / n.alleles
    return(freq.1)
  })
  
  # If requested, plot the site-frequency spectrum
  if(plot.sfs) {
    # Calculate minor allele frequency
    MAF <- apply(X = geno.mat, MARGIN = 2, FUN = function(snp) {
      # Find the total number of "1" alleles
      sum.1 <- sum(snp)
      # Find the total number of alleles
      n.alleles <- 2 * length(snp)
      # Calculate frequency
      freq.1 <- sum.1 / n.alleles
      maf <- min(freq.1, (1 - freq.1))
      return(maf)
    })
    
    # Plot
    hist(x = MAF,
         main = "Site Frequency Spectrum",
         xlab = "Minor Allele Frequency",
         ylab = "Count",
         xlim = c(0,0.5))
  }
  
  # Return the data
  return(freq.alleles)
  
} # Close the function

# # Define a function to filter a genotype matrix
# filter.genos <- function(geno.mat,
#                          max.marker.missing,
#                          max.marker.het,
#                          min.MAF,
#                          max.entry.het,
#                          max.entry.missing
# ) {
#   
#   # Deal with input
#   geno.mat <- as.matrix(geno.mat)
#   # Add 1 to get the 0, 1, 2 format
#   geno.mat.recode <- geno.mat + 1
#   
#   # Filter on markers first
#   # Calculate minor allele frequency
#   MAF <- apply(X = geno.mat.recode, MARGIN = 2, FUN = function(snp) {
#     freq <- sum(snp, na.rm = T) / (2 * length(snp))
#     maf <- min(freq, (1 - freq))
#     return(maf)
#   })
#   # Calculate marker heterozygosity
#   marker.het <- apply(X = geno.mat.recode, MARGIN = 2, FUN = function(snp) {
#     het.freq <- sum(snp == 1, na.rm = T) / length(snp)
#     return(het.freq)
#   })
#   # Calculate marker missingness
#   marker.missing <- apply(X = geno.mat.recode, MARGIN = 2, FUN = function(snp) {
#     missing.freq <- sum(is.na(snp)) / length(snp)
#     return(missing.freq)
#   })
#   
#   # Filter
#   to.keep <- which( MAF > min.MAF & marker.het < max.marker.het & marker.missing < max.marker.missing )
#   geno.mat.marker.filt <- geno.mat[,to.keep]
#   geno.mat.recode.marker.filt <- geno.mat.marker.filt + 1
#   
#   # Filter on entries next
#   # Calculate entry heterozygosity
#   entry.het <- apply(X = geno.mat.recode.marker.filt, MARGIN = 1, FUN = function(entry) {
#     het.freq <- sum(entry == 1, na.rm = T) / length(entry)
#     return(het.freq)
#   })
#   # Calculate entry missingness
#   entry.missing <- apply(X = geno.mat.recode.marker.filt, MARGIN = 1, FUN = function(entry) {
#     missing.freq <- sum(is.na(entry)) / length(entry)
#     return(missing.freq)
#   })
#   
#   # Filter
#   to.keep <- which( entry.het < max.entry.het & entry.missing < max.entry.missing )
#   geno.mat.filt <- geno.mat.marker.filt[to.keep,]
#   
#   # Return the matrix
#   return(geno.mat.filt)
# }

# Define a function for creating a genome
## This function is updated to allow for chromosomes to have unequal numbers of 
## loci / markers
make.genome <- function(n.chr, # An integer specifying the number of chromosomes
                        chr.len, # A numeric vector specifiying length of each chromosome.
                        n.chr.snps, # A numeric vector specifying the number of loci on each chromosome
                        genetic.map # A list of genetic map positions for loci on each chromosome
) {
  
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

# Define a function to define the trait architecture
trait.architecture <- function(genome,
                               n.QTL, # Number of desired QTL. If n.QTL / n.chr is non-integer, some QTL are given an effect of 0
                               qtl.index = NULL, # A list of QTL indices per chromosome. If NULL, indicies are randomly sampled
                               qtl.dom.index = NULL, # The index of QTL with dominance. If NULL, no QTL show dominance
                               qtl.perf.index = NULL, # The index of perfect loci (i.e. the SNP is the QTL). If NULL, there are no perfect markers
                               qtl.add.eff = "normal", # The additive effects of QTL. Can be "normal", "geometric", or a vector of same length as qtl.ids
                               qtl.dom.eff = NULL # The dominance effects of QTL. Can be "NULL" (no dominance) or a vector of same length as qtl.dom.ids
                               ) {
  
  # Deal with input
  # Find the number of chromosomes
  n.chr <- length(genome)
  
  # QTL.per.chr <- n.QTL / n.chr
  # # If the quotient is not an integer, round up
  # if (!is.integer(QTL.per.chr)) {
  #   QTL.per.chr <- ceiling(QTL.per.chr)
  # }
  # # Find the number of total QTL
  # n.tot.QTL = QTL.per.chr * n.chr
  
  # If the qtl.ids are NULL, randomly sample the loci
  if (is.null(qtl.index)) {
    
    # Create a list to determine the number of qtl for each chromosome. This will
    ## essentially make the distribution of qtl even across the chromosome indices
    qtl.per.chr <- sapply(X = split(x = 1:n.QTL, f = cut(1:n.QTL, breaks = n.chr)), FUN = length)
    
    # Create an empty list
    qtl.index.per.chr <- list()
    
    # Loop over the chromosomes in the genome
    for (p in 1:n.chr) {
      
      # Pull out the number of QTL designated for the chromosome
      n.qtl.p <- qtl.per.chr[p]
      
      # Pull out the number of loci on the chromsome and create an index
      loci.index.p <- seq(1, slot(genome[[p]], "num.snp.chr"))
      
      # Sample the index and add to the list
      qtl.index.per.chr[[p]] <- sort(sample(x = loci.index.p, size = n.qtl.p))
    }
    
    # If the list is specified, check it
  } else {
    
    if (length(qtl.index) != n.chr) stop("The length of the qtl.index list is not the same as the number of chromosomes.")
    # Rename it
    qtl.index.per.chr <- qtl.index
    # Figure out the QTL per chromosomes
    qtl.per.chr <- sapply(X = qtl.index.per.chr, FUN = length)
  }
  
  # QTL dominance IDs
  if (is.null(qtl.dom.index)) {
    # Create a list of NULLs
    qtl.dom.index.per.chr <- sapply(X = 1:n.chr, FUN = function(i) NULL)
    
  } else { # Otherwise check the list
    if (length(qtl.dom.index) != n.chr) stop("The length of the qtl.dom.index list is not the same as the number of chromsomes.")
    
    # Rename it
    qtl.dom.index.per.chr <- qtl.dom.index
  }
  
  # Perfect QTL IDs
  if (is.null(qtl.perf.index)) {
    # Create a list of NULLs
    qtl.perf.index.per.chr <- sapply(X = 1:n.chr, FUN = function(i) NULL)
    
  } else { # Otherwise check the list
    if (length(qtl.perf.index) != n.chr) stop("The length of the qtl.perf.index list is not the same as the number of chromosomes.")
    
    # Rename it
    qtl.perf.index.per.chr <- qtl.perf.index
  }
  
  # Assign qtl effects
  if (is.character(qtl.add.eff)) {
    
    # Draw from standard normal
    if (qtl.add.eff == "normal") {
      qtl.add.eff <- abs(rnorm(n = n.QTL, mean = 0, sd = 1))
      # Randomly assign negative values to the qtl effects - this corresponds to the value of the 1 allele
      qtl.add.eff <- qtl.add.eff * sample(c(1,-1), n.QTL, replace = T)
      
    } else { # Draw from geometric series
      
      if (qtl.add.eff == "geometric") {
        a = (n.QTL - 1) / (n.QTL + 1)
        qtl.add.eff <- sample(a^(1:n.QTL))
        # Randomly assign negative values to the qtl effects - this corresponds to the value of the 1 allele
        qtl.add.eff <- qtl.add.eff * sample(c(1,-1), n.QTL, replace = T)
        
        # Break up the effects into chromosomes with the same number of elements
        ## as QTL on those chromosomes
        qtl.add.eff.per.chr <- split(qtl.add.eff, rep(1:n.chr, qtl.per.chr))
        
      } else {
        # Otherwise, make sure the length of the additive effects is the same as the number of chr
        if (length(qtl.add.eff) != n.chr) stop("The length of the QTL effects is not the same as the number of chromosomes")
      }}}
  
  # Dominance effects
  if (is.null(qtl.dom.eff)) {
    # Create a list of NULLs
    qtl.dom.eff.per.chr <- sapply(X = 1:n.chr, FUN = function(i) NULL)
    
  } else { # Otherwise check it
    if (length(qtl.dom.eff) != n.chr) stop("The length of the qtl.dom.eff list is not the same as the number of chromosomes.")
    
    # Make sure each element in the qtl.dom.eff list is the same as the qtl.dom.index.per.chr list
    dom.elements.same <- sapply(X = 1:n.chr, FUN = function(i) length(qtl.dom.eff[[i]]) == length(qtl.dom.index.per.chr[[i]]) )
    
    if (!all(dom.elements.same)) stop("One or more of the elements in the qtl.dom.eff list are not the same length as the corresponding qtl.dom.index list.")
  }


  
  # Apply a function over each chromosome in the genome
  genome <- sapply(X = 1:n.chr, FUN = function(i) {
    
    # Add to the ith position in the genome list a revised chromosome genome with
    ## the qtl positions and effects
    genome[[i]] <- hypredNewQTL(genome[[i]],
                                new.id.add = qtl.index.per.chr[[i]],
                                new.id.dom = qtl.dom.index.per.chr[[i]],
                                new.id.per.mar = qtl.perf.index.per.chr[[i]],
                                new.eff.add = qtl.add.eff.per.chr[[i]],
                                new.eff.dom = qtl.dom.eff.per.chr[[i]] )
  })
  
  # Return the new genome
  return(genome)
  
} # Close the function



# Define a function to simulate recombination. This function will use the
## hypredRecombine function as a base, but will apply that function over
## the chromosomes in the genome list
recombine <- function(genome,
                      haploid.genomeA,
                      haploid.genomeB,
                      mutate = TRUE,
                      mutation.rate.snp,
                      mutation.rate.qtl,
                      block = FALSE
){
  
  # Pull out the number of chromosomes
  n.chr <- length(genome)
  # Pull out the number of loci per chromsosome
  loci.per.chr <- sapply(X = genome, FUN = function(chromosome) length(slot(chromosome, "pos.snp")))
  
  # Split each haploid genome into chromosomes
  haploid.genomeA.split <- split(haploid.genomeA, rep(1:n.chr, loci.per.chr))
  haploid.genomeB.split <- split(haploid.genomeB, rep(1:n.chr, loci.per.chr))
  
  # Apply a function over each chromosome
  gamete.list <- sapply(X = 1:n.chr, FUN = function(i) {
    
    hypredRecombine(genome[[i]],
                    genomeA = haploid.genomeA.split[[i]],
                    genomeB = haploid.genomeB.split[[i]],
                    mutate = mutate,
                    mutation.rate.snp = mutation.rate.snp,
                    mutation.rate.qtl = mutation.rate.qtl,
                    block = block) })
  
  # Concatenate the gamete
  return(do.call("c", gamete.list))
  
} # Close the function
  

# Define a function to optimize the training set based on the
## CDmean maximization algorithm. This code is taken from Rincent et al 2012
TP.optimization <- function(genos, # genotype matrix of all candidate individuals
                            phenotyped.lines, # Character of line name to be considered "phenotyped"
                            unphenotyped.lines, # Character of line name to be considered "unphenotyped"
                            n.TP, # The size of the desired TP
                            V_e, # V_e for calculating lambda
                            V_a, # V_a for calculating lambda
                            optimization = c("PEVmean", "CDmean"),
                            max.iter = 800,
                            use.subset = FALSE # Logical indicating whether a subset of the whole should be used for measuring the contrats. This saves time but may be innacurate
                            ) {
  
  # Define a function to create a matrix of contrasts between the individuals not
  ## in the training population (the candidates) and the mean of the
  ## whole population. This is taken from Rincent et al 2012
  make.contrast <- function(unphenotyped.index, # Index of the unphenotyped individuals in the A.mat
                                  n.total) {
    # Number of unphenotyped individuals
    n.unphenotyped <- length(unphenotyped.index)
    
    # The positive value of the contrast
    pos.contrast = 1 - (1 / n.total)
    # The negative value of the contrast
    neg.contrast = -1 / n.total
    
    # Make the total constrast matrix
    contrast.mat <- matrix(neg.contrast, n.total, n.unphenotyped)
    # Fill in the positive values
    for (i in 1:n.unphenotyped) {
      contrast.mat[unphenotyped.index[i], i] <- pos.contrast
    }
    
    return(contrast.mat)
  }

  # Calculate the relationship matrix
  # The first set of lines is the phenotyped and the second set is the unphenotyped
  A <- A.mat(genos, 0, 1)
  
  # Find the index of the phenotyped and unphenotyped lines in the A matrix
  phenotyped.index <- which(row.names(A) %in% phenotyped.lines)
  unphenotyped.index <- which(row.names(A) %in% unphenotyped.lines)
  
  # Set the total
  # If using a subset, the total is size of the training set + the number of unphenotyped lines
  if (use.subset) {
    n.total = n.TP + length(unphenotyped.index)
  } else {
    # Otherwise it is the number of lines in the A matrix
    n.total = nrow(A)
  }
  
  # Calculate lambda using the provided variances
  lambda = V_e / V_a
  
  ## Begin the algorithm
  # The algorithm is the same for both PEVmean and CDmean (exchange), so we have two paths
  ## depending on which is desired
  if (optimization == "PEVmean") {
    
    # Set the CDmean to NA
    CDmean.save = NA
    CD.vector = NA
    
    ### Initialize the optimization algorithm
    ## Separate paths based on whether a subset will be used
    if (use.subset) {
      
      # Create design matricies
      I <- diag(n.TP) # Identity matrix of TP size
      X <- matrix(1, n.TP, 1) # X matrix
      M <- I - ( X %*% solve(t(X) %*% X) %*% t(X) )
      
      # Randomly sample the "phenotyped" lines to start the TP
      ## These will serve as the first sample of phenotyped lines
      start.OTP <- sort(sample(phenotyped.index, n.TP))
      # Save the sample
      save.OTP <- start.OTP
    
      # Find the "phenotyped" lines that were not included, but could be included
      candidates.OTP <- setdiff(phenotyped.index, start.OTP)
  
      # When using a subset (to save time)
      Z <- diag(n.total)[1:length(start.OTP),] # This matrix will remain constant
      # Create a vector of the indicies of the sampled training set and the unphenotyped lines
      subset.index <- c(start.OTP, unphenotyped.index)
      # Subset the A.matrix
      A1 <- A[subset.index, subset.index]
      
      # The index of the unphenotyped lines becomes the last n rows of the subsetted A matrix, where n
      ## is the number of unphenotyped lines
      unphenotyped.index.A1 <- setdiff(1:nrow(A1), 1:length(start.OTP))
      
      # Try to find the inverse of candidate.A.mat, if not add ridge
      A1.inv <- try(solve(A1), silent = T)
      # Assess whether it was an error
      if (class(A1.inv) == "try-error") {
        r = 1e-6
        ridge <- diag(nrow(A1)) * r
        A1 = A1 + ridge
        # Recalculate the inverse
        A1.inv <- solve(A1)
      }
      
      # Contrasts matrix
      ## The dimenstions should be row = n.total and column = n.unphenotyped (aka the parents)
      c.mat <- make.contrast(unphenotyped.index = unphenotyped.index.A1, n.total = n.total) # This matrix will also remain constant
      
      # Calculate the inital PEVmean of the set
      inv.C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
      numerator = ( t(c.mat) %*% inv.C_22 %*% c.mat )
      denominator = ( t(c.mat) %*% c.mat )
      PEV.mat = numerator / denominator
      
      PEV = diag(PEV.mat) * V_e
      # Calculate the PEVmean
      PEVmean.save = mean(PEV)
      
      # Create an empty vector to store PEVmean values
      PEV.vector <- numeric()
      
      iter.counter = 1 # Iteration counter
      # While loop
      while(iter.counter <= max.iter) {
        
        # Add one to the count
        iter.counter = iter.counter + 1
        
        # Randomly remove one individual from the training set
        train.to.remove <- sample(start.OTP, 1)
        # Randomly choose one individual from the candidate set to add
        candidate.to.add <- sample(candidates.OTP, 1)
        
        # Create the new training set
        new.OTP <- sort(c(candidate.to.add, start.OTP[start.OTP != train.to.remove]))
        
        # Create the new candidate set
        candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        
        # Create a vector of the indicies of the sampled training set and the unphenotyped lines
        subset.index <- c(new.OTP, unphenotyped.index)
        # Subset the A.matrix
        A1 <- A[subset.index, subset.index]
        
        # Try to find the inverse of candidate.A.mat, if not add ridge
        A1.inv <- try(solve(A1), silent = T)
        # Assess whether it was an error
        if (class(A1.inv) == "try-error") {
          r = 1e-6
          ridge <- diag(nrow(A1)) * r
          A1 = A1 + ridge
          # Recalculate the inverse
          A1.inv <- solve(A1)
        }
        
        inv.C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
        numerator = ( t(c.mat) %*% inv.C_22 %*% c.mat )
        denominator = ( t(c.mat) %*% c.mat )
        PEV.mat = numerator / denominator
        
        PEV = diag(PEV.mat) * V_e
        PEVmean.new <- mean(PEV)
        
        # If the new PEVmean is lower, accept, if not reject
        if (PEVmean.new < PEVmean.save) {
          start.OTP <- new.OTP
          PEVmean.save <- PEVmean.new
        }
        
        # Record the PEVmean of the new PEVmean (if accepted) or the previous (if not)
        PEV.vector[iter.counter - 1] <- PEVmean.save
        
      } # Close the loop
      
    } else { # Path for using the whole set
      
      # Create design matricies
      I <- diag(n.TP) # Identity matrix of TP size
      X <- matrix(1, n.TP, 1) # X matrix
      M <- I - ( X %*% solve(t(X) %*% X) %*% t(X) )
      
      # Randomly sample the "phenotyped" lines to start the TP
      ## These will serve as the first sample of phenotyped lines
      start.OTP <- sort(sample(phenotyped.index, n.TP))
      # Save the sample
      save.OTP <- start.OTP
      
      # Find the "phenotyped" lines that were not included, but could be included
      candidates.OTP <- setdiff(phenotyped.index, start.OTP)
      
      # When using the whole candidates
      Z <- matrix(0, n.TP, n.total)
      # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
      for (i in 1:length(start.OTP)) { Z[i, start.OTP[i] ] = 1 } 
      # Use the whole A.mat
      A1 <- A

      # Try to find the inverse of candidate.A.mat, if not add ridge
      A1.inv <- try(solve(A1), silent = T)
      # Assess whether it was an error
      if (class(A1.inv) == "try-error") {
        r = 1e-6
        ridge <- diag(nrow(A1)) * r
        A1 = A1 + ridge
        # Recalculate the inverse
        A1.inv <- solve(A1)
      }
      
      # Contrasts matrix
      ## The dimenstions should be row = n.total and column = n.unphenotyped (aka the parents)
      c.mat <- make.contrast(unphenotyped.index = unphenotyped.index, n.total = n.total)

      # Calculate the inital PEVmean of the set
      inv.C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
      numerator = ( t(c.mat) %*% inv.C_22 %*% c.mat )
      denominator = ( t(c.mat) %*% c.mat )
      PEV.mat = numerator / denominator
      
      PEV = diag(PEV.mat) * V_e
      # Calculate the PEVmean
      PEVmean.save = mean(PEV)
      
      # Create an empty vector to store PEVmean values
      PEV.vector <- numeric()
      
      iter.counter = 1 # Iteration counter
      # While loop
      while(iter.counter <= max.iter) {
        
        # Add one to the count
        iter.counter = iter.counter + 1
        
        # Randomly remove one individual from the training set
        train.to.remove <- sample(start.OTP, 1)
        # Randomly choose one individual from the candidate set to add
        candidate.to.add <- sample(candidates.OTP, 1)
        
        # Create the new training set
        new.OTP <- sort(c(candidate.to.add, start.OTP[start.OTP != train.to.remove]))
        
        # Create the new candidate set
        candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        # When using the whole candidates
        Z <- matrix(0, n.TP, n.total)
        # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
        for (i in 1:length(new.OTP)) { Z[i, new.OTP[i] ] = 1 } 

        # Calculate the inital PEVmean of the set
        # No need to remake the contrast matrix because the unphenotyped lines remain unchanged
        inv.C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
        numerator = ( t(c.mat) %*% inv.C_22 %*% c.mat )
        denominator = ( t(c.mat) %*% c.mat )
        PEV.mat = numerator / denominator
        
        PEV = diag(PEV.mat) * V_e
        PEVmean.new <- mean(PEV)
        
        # If the new PEVmean is lower, accept, if not reject
        if (PEVmean.new < PEVmean.save) {
          start.OTP <- new.OTP
          PEVmean.save <- PEVmean.new
        }
        
        # Record the PEVmean of the new PEVmean (if accepted) or the previous (if not)
        PEV.vector[iter.counter - 1] <- PEVmean.save
        
      } # Close the loop
      
    } # Close the use.subset if statement
    
  } # Close the PEVmean optimization if statement
  
  # New if statement for CDmean
  if(optimization == "CDmean") {
    
    # Set the PEVmean to NA
    PEVmean.save = NA
    PEV.vector = NA
    
    ### Initialize the optimization algorithm
    ## Separate paths based on whether a subset will be used
    if (use.subset) {
      
      # Create design matricies
      I <- diag(n.TP) # Identity matrix of TP size
      X <- matrix(1, n.TP, 1) # X matrix
      M <- I - ( X %*% solve(t(X) %*% X) %*% t(X) )
      
      # Randomly sample the "phenotyped" lines to start the TP
      ## These will serve as the first sample of phenotyped lines
      start.OTP <- sort(sample(phenotyped.index, n.TP))
      # Save the sample
      save.OTP <- start.OTP
      
      # Find the "phenotyped" lines that were not included, but could be included
      candidates.OTP <- setdiff(phenotyped.index, start.OTP)
      
      # When using a subset (to save time)
      Z <- diag(n.total)[1:length(start.OTP),] # This matrix will remain constant
      # Create a vector of the indicies of the sampled training set and the unphenotyped lines
      subset.index <- c(start.OTP, unphenotyped.index)
      # Subset the A.matrix
      A1 <- A[subset.index, subset.index]
      
      # The index of the unphenotyped lines becomes the last n rows of the subsetted A matrix, where n
      ## is the number of unphenotyped lines
      unphenotyped.index.A1 <- setdiff(1:nrow(A1), 1:length(start.OTP))
      
      # Try to find the inverse of candidate.A.mat, if not add ridge
      A1.inv <- try(solve(A1), silent = T)
      # Assess whether it was an error
      if (class(A1.inv) == "try-error") {
        r = 1e-6
        ridge <- diag(nrow(A1)) * r
        A1 = A1 + ridge
        # Recalculate the inverse
        A1.inv <- solve(A1)
      }
      
      # Contrasts matrix
      ## The dimenstions should be row = n.total and column = n.unphenotyped (aka the parents)
      c.mat <- make.contrast(unphenotyped.index = unphenotyped.index.A1, n.total = n.total) # This matrix will also remain constant
      
      # Calculate the inital CDmean of the set
      inv.C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
      numerator = ( t(c.mat) %*% (A1 - (lambda * inv.C_22)) %*% c.mat )
      denominator = t(c.mat) %*% A1 %*% c.mat
      CD.mat = numerator / denominator
      
      CD = diag(CD.mat)
      # Calculate the CDmean
      CDmean.save = mean(CD)
      
      # Create an empty vector to store PEVmean values
      CD.vector <- numeric()
      
      iter.counter = 1 # Iteration counter
      # While loop
      while(iter.counter <= max.iter) {
        
        # Add one to the count
        iter.counter = iter.counter + 1
        
        # Randomly remove one individual from the training set
        train.to.remove <- sample(start.OTP, 1)
        # Randomly choose one individual from the candidate set to add
        candidate.to.add <- sample(candidates.OTP, 1)
        
        # Create the new training set
        new.OTP <- sort(c(candidate.to.add, start.OTP[start.OTP != train.to.remove]))
        
        # Create the new candidate set
        candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        
        # Create a vector of the indicies of the sampled training set and the unphenotyped lines
        subset.index <- c(new.OTP, unphenotyped.index)
        # Subset the A.matrix
        A1 <- A[subset.index, subset.index]
        
        # Try to find the inverse of candidate.A.mat, if not add ridge
        A1.inv <- try(solve(A1), silent = T)
        # Assess whether it was an error
        if (class(A1.inv) == "try-error") {
          r = 1e-6
          ridge <- diag(nrow(A1)) * r
          A1 = A1 + ridge
          # Recalculate the inverse
          A1.inv <- solve(A1)
        }
        
        inv.C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
        numerator = ( t(c.mat) %*% (A1 - (lambda * inv.C_22)) %*% c.mat )
        denominator = t(c.mat) %*% A1 %*% c.mat
        CD.mat = numerator / denominator
        
        CD = diag(CD.mat)
        # Find the new mean
        CDmean.new <- mean(CD)
        
        # If the new CDmean is lower, accept, if not reject
        if (CDmean.new > CDmean.save) {
          start.OTP <- new.OTP
          CDmean.save <- CDmean.new
        }
        
        # Record the CDmean of the new CDmean (if accepted) or the previous (if not)
        CD.vector[iter.counter - 1] <- CDmean.save
        
      } # Close the loop
      
    } else { # Path for using the whole set
      
      # Create design matricies
      I <- diag(n.TP) # Identity matrix of TP size
      X <- matrix(1, n.TP, 1) # X matrix
      M <- I - ( X %*% solve(t(X) %*% X) %*% t(X) )
      
      # Randomly sample the "phenotyped" lines to start the TP
      ## These will serve as the first sample of phenotyped lines
      start.OTP <- sort(sample(phenotyped.index, n.TP))
      # Save the sample
      save.OTP <- start.OTP
      
      # Find the "phenotyped" lines that were not included, but could be included
      candidates.OTP <- setdiff(phenotyped.index, start.OTP)
      
      # When using the whole candidates
      Z <- matrix(0, n.TP, n.total)
      # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
      for (i in 1:length(start.OTP)) { Z[i, start.OTP[i] ] = 1 } 
      # Use the whole A.mat
      A1 <- A
      
      # Try to find the inverse of candidate.A.mat, if not add ridge
      A1.inv <- try(solve(A1), silent = T)
      # Assess whether it was an error
      if (class(A1.inv) == "try-error") {
        r = 1e-6
        ridge <- diag(nrow(A1)) * r
        A1 = A1 + ridge
        # Recalculate the inverse
        A1.inv <- solve(A1)
      }
      
      # Contrasts matrix
      ## The dimenstions should be row = n.total and column = n.unphenotyped (aka the parents)
      c.mat <- make.contrast(unphenotyped.index = unphenotyped.index, n.total = n.total)
      
      # Calculate the inital PEVmean of the set
      inv.C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
      numerator = ( t(c.mat) %*% (A1 - (lambda * inv.C_22)) %*% c.mat )
      denominator = t(c.mat) %*% A1 %*% c.mat
      CD.mat = numerator / denominator
      
      CD = diag(CD.mat)
      # Calculate the CDmean
      CDmean.save = mean(CD)
      
      # Create an empty vector to store CDmean values
      CD.vector <- numeric()
      
      iter.counter = 1 # Iteration counter
      # While loop
      while(iter.counter <= max.iter) {
        
        # Add one to the count
        iter.counter = iter.counter + 1
        
        # Randomly remove one individual from the training set
        train.to.remove <- sample(start.OTP, 1)
        # Randomly choose one individual from the candidate set to add
        candidate.to.add <- sample(candidates.OTP, 1)
        
        # Create the new training set
        new.OTP <- sort(c(candidate.to.add, start.OTP[start.OTP != train.to.remove]))
        
        # Create the new candidate set
        candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        # When using the whole candidates
        Z <- matrix(0, n.TP, n.total)
        # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
        for (i in 1:length(new.OTP)) { Z[i, new.OTP[i] ] = 1 } 
        
        # Calculate the inital PEVmean of the set
        # No need to remake the contrast matrix because the unphenotyped lines remain unchanged
        inv.C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
        numerator = ( t(c.mat) %*% (A1 - (lambda * inv.C_22)) %*% c.mat )
        denominator = t(c.mat) %*% A1 %*% c.mat
        CD.mat = numerator / denominator
        
        CD = diag(CD.mat)
        # Find the new mean
        CDmean.new <- mean(CD)
        
        # If the new CDmean is lower, accept, if not reject
        if (CDmean.new > CDmean.save) {
          start.OTP <- new.OTP
          CDmean.save <- CDmean.new
        }
        
        # Record the CDmean of the new CDmean (if accepted) or the previous (if not)
        CD.vector[iter.counter - 1] <- CDmean.save
        
      } # Close the loop
      
    } # Close the use.subset if statement
    
  } # Close the CD optimization if statement

  ### Output the optimized TP
  # Extract line names from A.mat
  line.names <- row.names(A)
  
  OTP.lines <- line.names[start.OTP]
  
  output.list <- list(optimized.lines = OTP.lines,
                      PEVmean = list(PEVmean = PEVmean.save,
                                     PEV.vector = PEV.vector),
                      CDmean = list(CDmean = CDmean.save,
                                    CD.vector = CD.vector) )
                                    
  
  return(output.list)
} # Close the function


# Define a function to take a gamete matrix and add additional columns to reflect qtl
augment.gamete.mat <- function(genome, # The object of class "hypredGenome"
                               gamete.mat # A 2n x m matrix of gamete loci values for each line
                               ) {
  
  # Pull out the qtl and snp loci from the genome
  qtl.ind <- slot(genome, "pos.add.qtl")$ID
  snp.ind <- setdiff( 1:length(slot(genome, "pos.snp")), qtl.ind )
  
  # Pull out the marker names
  marker.names <- colnames(gamete.mat)
  # Create QTL names
  qtl.names <- paste("QTL", 1:length(qtl.ind), sep = ".")
  
  # Find the total number of sites
  tot.loci <- length(c(snp.ind, qtl.ind))
  
  # Make sure that the number of snps and qtl doesn't already equal the number of loci in the gamete.mat
  if (tot.loci == ncol(gamete.mat)) {
    # Issue a warning
    warning("The total number of loci (SNPs and QTL) in the genome already equals the number of loci in the gamete matrix. 
            The function will return the gamete matrix unchanged.")
    # Return the gamete matrix unchanged
    return(gamete.mat)
  }
  
  # Find the number of haploid genomes
  n.gametes <- nrow(gamete.mat)
  # Make a matrix of the appropriate size
  new.gamete.mat <- matrix(NA, nrow = n.gametes, ncol = tot.loci)
  
  # Add the SNP genotypes to the snp posistions in the new gamete matrix
  new.gamete.mat[,snp.ind] <- gamete.mat
  
  # Randomly generate allele states for the qtl 
  qtl.allele.states <- replicate(length(qtl.ind), rep(round(runif(n.gametes / 2, 0, 1)), each = 2) )
  new.gamete.mat[,qtl.ind] <- qtl.allele.states
  
  # Add dimnames
  row.names(new.gamete.mat) <- row.names(gamete.mat)
  
  colnames.vec <- character()
  colnames.vec[snp.ind] <- marker.names
  colnames.vec[qtl.ind] <- qtl.names
  colnames(new.gamete.mat) <- colnames.vec
  
  # Return the matrix
  return(new.gamete.mat)
  
} # Close the function


# Define a function to measure the LD between adjacent markers and QTL on the same
## chromosome. The function will measure the LD between all polymorphic QTL and markers by finding
## the correlation between that QTL and the markerrs
measure.LD <- function(genome, # Genome object
                       genos # n x m matrix of genotype data with QTL
                       ) {
  
  # Pull out the number of snps per chromosome
  loci.per.chr <- c(0, sapply(genome, function(chromosome) length(slot(chromosome, "pos.snp"))))
  
  # Make empty lists
  pos.qtl <- list()

  # Find the index of the qtl in the combined matrix  
  for (i in 1:length(genome)) {
    
    # QTL index for the ith chromosome
    qtl.index.i <- slot(genome[[i]], "pos.add.qtl")$ID
    # Adjust for the previous chromosomes and add to the vector
    pos.qtl[[i]] <- sum(loci.per.chr[1:i]) + qtl.index.i
  }
    
  # Collapse to vector
  pos.qtl <- do.call("c", pos.qtl)
  # Figure the positions of the SNPs
  pos.snps <- setdiff(1:ncol(genos), pos.qtl)
  
  # Determine the index of the polymorphic qtl
  poly.qtl.index <- pos.qtl[apply(X = genos[,pos.qtl], MARGIN = 2, FUN = function(locus) length(unique(locus)) > 1)]
  # Determine the index of the polymorphic markers (monomorphic markers will return a NA for the correlation)
  poly.snp.index <- pos.snps[apply(X = genos[,pos.snps], MARGIN = 2, FUN = function(locus) length(unique(locus)) > 1)]
  
  # Apply the function over the genotype matrix of just the polymorphic qtl
  poly.qtl.marker.LD <- sapply(X = poly.qtl.index, FUN = function(qtl.id) {
    
    # Subset the geno matrix
    qtl <- genos[,qtl.id]
    
    # Measure LD as the correlation between the qtl genotypes and the marker genotypes
    apply(X = genos[,poly.snp.index], MARGIN = 2, FUN = function(snp)
      # Correlate
      cor(qtl, snp) ) })
  
  # Return the results
  return(t(poly.qtl.marker.LD))
  
} # Close the function
  
  # # Conver cM to M
  # sliding.window.M <- sliding.window.cM / 100
  # 
  # # Deal with input
  # if (!is.list(gametes)) {
  #   gamete.mat <- as.matrix(gametes)
  # } else {
  #   gamete.mat <- as.matrix(do.call("rbind", gametes))
  # }
  
  # # Pull out genome info
  # n.chr <- slot(genome, "num.chr")
  # 
  # # Pull out QTL and marker locations
  # loci.pos.M <- slot(genome, "pos.snp")
  # names(loci.pos.M) <- colnames(genos)
  
  # # Split loci positions into chromosomes
  # loci.pos.M.split <- split(x = loci.pos.M, cut(x = 1:length(loci.pos.M), breaks = n.chr))
  # 
  # # Define functions to calculate LD
  # LD <- function(g1, g2) {
  #   # Calculate frequency of 1 allele in the snp and qtl
  #   p.A <- sum(g1 == 1) / length(g1)
  #   p.a = 1 - p.A
  #   p.B <- sum(g2 == 1) / length(g2)
  #   p.b = 1 - p.B
  #   p.AB <- sum(apply(X = cbind(g1, g2), MARGIN = 1, FUN = function(allele.pair) all(allele.pair == c(1,1)) )) / length(g1)
  # 
  #   # Calculate LD
  #   # D
  #   D = p.AB - (p.A * p.B)
  #   
  #   # D'
  #   if (D == 0) {
  #     # If D is 0, so is D.prime
  #     D.prime = 0
  #   } else {
  #     if (D > 0) {
  #       D.max = min( (p.A * p.b), (p.a * p.B) )
  #     } else {
  #       D.max = max( -(p.A * p.B), -(p.a * p.b) )
  #     }
  #     # Calculate D prime
  #     D.prime = D / D.max
  #   }
  #   
  #   # r
  #   r = -D / sqrt( (p.A * p.a * p.B * p.b) )
  #   
  #   # Return a vector
  #   out.vec <- c(D, D.prime, r); names(out.vec) <- c("D", "D.prime", "r")
  #   return(out.vec)
  # }
  
  # # Define a function to split a vector of genetic positions by a certain window
  # split.pos <- function(pos.vector, window) {
  #   # Assign groups
  #   f <- as.factor(sapply(X = pos.vector, FUN = function(i) findInterval(x = c(i, i + window), vec = pos.vector))[2,])
  #   # Split
  #   return(split(x = pos.vector, f = f))
  # }
  #   
  # # Apply a function over the split loci list
  # pairwise.LD <- lapply(X = loci.pos.M.split, FUN = function(chr) {
  #   # Find the loci that are separated in distance by no more than the sliding window length
  #   loci.split <- split.pos(pos.vector = chr, window = sliding.window.M)
  #   # Remove groups of less than 2
  #   loci.split <- loci.split[sapply(X = loci.split, FUN = length) > 1]
  #   # Apply LD function over the groups
  #   do.call("cbind", sapply(X = loci.split, FUN = function(group) {
  #     # Generate combinations
  #     combn.mat <- data.frame(t(combn(x = names(group), 2)))
  #     # Apply over combinations
  #     combn.LD <- apply(X = combn.mat, MARGIN = 1, FUN = function(pair) {
  #       LD(g1 = gamete.mat[,pair[1]], g2 = gamete.mat[,pair[2]]) })
  #     # Add columns
  #     colnames(combn.LD) <- apply(X = combn.mat, MARGIN = 1, FUN = paste, collapse = ".")
  #     return(combn.LD)
  #   })) })
  # 
  # names(pairwise.LD) <- 1:n.chr
  # 
  # # Return list
  # return(list(pairwise.LD = pairwise.LD, sliding.window.cM = sliding.window.cM))
  



# Define a function to determine the mean relationship between TP individuals
## and selection candidates
measure.relationship <- function(TP.genos, # An n x m matrix of marker genotypes for the TP
                                 candidates.genos # A n x m matrix of marker genotypes for the selection candidates
                                 ) {
  
  # Pull out the names from each matrix
  TP.line.names <- row.names(TP.genos)
  candidate.line.names <- row.names(candidates.genos)
  
  # Combine the matricies
  all.genos <- rbind(TP.genos, candidates.genos)
  
  # Find the additive relationship matrix
  G <- A.mat(all.genos, min.MAF = 0, max.missing = 1)
  
  # Subset the G matrix for TP (rows) and candidates (columns)
  G.sub <- G[TP.line.names, candidate.line.names]
  
  # Calculate the mean relationship of the whole TP to every candidate
  mu.rel.candidate <- colMeans(G.sub)
  # Find the mean of the mean relationship. This is essentially the mean relationship
  ## of the whole TP to all selection candidates
  mu.rel <- mean(mu.rel.candidate)
  
  # Return the mean
  return(mu.rel)

} # Close the function
  
  
      
# # To estimate variance, we will use a generalized linear mixed model
# # Genos and gxe will be fit as random effects, and environment will be fit as a fixed effect
# # First develop factors
# g_i = factor(rep(1:n.geno, n.env * n.rep))
# e_j = factor(rep(1:n.env, each = n.geno, n.rep))
# r_k = factor(rep(1:n.rep, each = n.geno * n.env))
# 
# # If no model is called, don't estimate variances
# if (is.null(model)) {
#   V_g.hat <- NA
#   V_residual.hat <- NA
#   h2.hat <- NA
#   V_p.hat <- NA
#   
# } else {
#   
#   # If the model is anova, run an anova
#   if (model == "anova") {
#     
#     p <- as.vector(y)
#     
#     # Phenotypic variance
#     V_p.hat <- var(p)
#     
#     # With more than one rep and more than one environment
#     if (all(n.rep > 1, n.env > 1)) {
#       p.anova <- anova(aov(p ~ g_i + e_j + r_k %in% e_j + g_i:e_j))
#       
#       # Find the mean squares
#       MS_g <- p.anova['g_i',]$`Mean Sq`
#       MS_ge <- p.anova['g_i:e_j',]$'Mean Sq'
#       MS_error <- p.anova['Residuals',]$`Mean Sq`
#       
#       # Calculate variance
#       V_g.hat <- (MS_g - MS_ge) / (n.env * n.rep)
#       V_ge.hat <- (MS_ge - MS_error) / n.rep
#       V_residual.hat <- MS_error
#       
#       # Calculate heritability
#       h2.hat = V_g.hat / (V_g.hat + (V_ge.hat / n.env) + (V_residual.hat / (n.env * n.rep)))
#       
#     } else {
#       # With 1 rep
#       if (n.rep == 1) {
#         p.anova <- anova(aov(p ~ g_i + e_j))
#         
#         # Find the mean squares
#         MS_g <- p.anova['g_i',]$`Mean Sq`
#         MS_error <- p.anova['Residuals',]$`Mean Sq`
#         
#         # Calculate variances
#         V_g.hat = (MS_g - MS_error) / (n.env * n.rep)
#         V_ge.hat <- NA
#         V_residual.hat = MS_error
#         
#         # Calculate heritability
#         h2.hat = V_g.hat / (V_g.hat + (V_residual.hat / (n.env * n.rep)))
#       }
#       
#       # With 1 env
#       if (n.env == 1) {
#         p.anova <- anova(aov(p ~ g_i + r_k))
#         
#         # Find the mean squares
#         MS_g <- p.anova['g_i',]$`Mean Sq`
#         MS_error <- p.anova['Residuals',]$`Mean Sq`
#         
#         # Calculate variance
#         V_g.hat = (MS_g - MS_error) / (n.env * n.rep)
#         V_ge.hat <- NA
#         V_residual.hat = MS_error
#         
#         # Calculate heritability
#         h2.hat = V_g.hat / (V_g.hat + (V_residual.hat / (n.env * n.rep)))
#       }
#     }
#   }
#   
#   # Use the mixed model if called on
#   if (model == "mixed") {
#     
#     # Stack the data 
#     y.stack <- data.frame(value = stack(as.data.frame(y))[,1], env = e_j, geno = g_i, rep = r_k)
#     
#     if (all(n.rep > 1, n.env > 1)) {
#       stop("This option is not yet developed")
#       # Solve a generalized linear mixed model
#       glm.fit <- lmer(value ~ -1 + env + (1|geno) + (geno|env), data = y.stack)
#       
#       # Extract variance estimates
#       V_g.hat <- VarCorr(glm.fit)$geno[1]
#       V_e.hat <- var(fixef(glm.fit))
#       V_residual.hat <- attr(VarCorr(glm.fit), "sc") ^ 2
#       
#     } else {
#       # With 1 rep
#       if (n.rep == 1) {
#         # Solve a generalized linear mixed model
#         glm.fit <- lmer(value ~ -1 + (1|geno) + env, data = y.stack)
#         
#         # Extract variance estimates
#         V_g.hat <- VarCorr(glm.fit)$geno[1]
#         V_e.hat <- var(fixef(glm.fit))
#         V_residual.hat <- attr(VarCorr(glm.fit), "sc") ^ 2
#         
#         # Estimate heritability
#         h2.hat = V_g.hat / (V_g.hat + (V_residual.hat / n.env))
#       }
#       
#       if (n.env == 1) {
#         glm.fit <- lmer(value ~ -1 + (1|geno), data = y.stack)
#         
#         # Extract variance estimates
#         V_g.hat <- VarCorr(glm.fit)$geno[1]
#         V_e.hat <- NA
#         V_residual.hat <- attr(VarCorr(glm.fit), "sc") ^ 2
#         
#         # Estimate heritability
#         h2.hat <- V_g.hat / (V_g.hat + (V_residual.hat / n.rep))
#       }
#     }
#   }
# } # Close the model null else statement



  
  
  
  
  
  
  
  
  
  