## Using "hypred" package for simulations

library(hypred)
library(EMMREML)
library(rrBLUP)

setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Side Projects/Simulations/Hypred simulations/")

## Demonstration according to the vignette
# 1. 250 F2 individuals from two founder lines are created
# 2. this population is then random mated for 50 generations
# 3. after this, the best individual is selected according to its phenotype (h2 =
#                                                                              0:5)
# 4. a DH population is then generated from this selected individual
# 5. this DH population is the training population for some genomic selection
# algorithm, hence the nal step is to generate the phenotypes and a design
# matrix

# First simulate the genome
# 2 chromosomes, each of 1.0 Morgans
# 100 SNP loci per chromosome
genome <- hypredGenome(num.chr = 2, len.chr = c(1.0, 1.0), num.snp.chr = 100)

# Assign a new genetic map where loci are evenly spaced
map <- rep(seq(0, 0.99, by = 0.01), 2)
genome <- hypredNewMap(genome, new.map = map)

# Return the genetic distance and expected recombination between sNPs
hypredSNPdist(genome, chromosome = 1, SNP1 = 1, SNP2 = 3)

# Assign some QTL
qtl.ids <- c(1, 20, 40, 60, 80, # chromosome 1
             101, 120, 140, 160, 180) # chromosome 2

# The third QTL in each chromosome has dominance effects
qtl.dom.ids <- c(40, # chr 1
                 140) # chr 2

# The second QTL in each chromosome is the perfect marker (i.e. the causative SNP)
per.mar.ids <- c(20, # chr 1
                 120) # chr 2

# Add QTL to the genome
genome <- hypredNewQTL(genome,
                       new.id.add = qtl.ids,
                       new.id.dom = qtl.dom.ids,
                       new.id.per.mar = per.mar.ids,
                       new.eff.add = rep(1, 10), # The additive effects are assumed to be 1
                       new.eff.dom = c(0.5, 0.5) # The dominance effects are assumed to be 0.5
                       )

# Summary
summary(genome)

# F2 base population
# Create the founder lines
founder1 <- hypredFounder(genome, prob.snp = 1)[1,]
founder2 <- hypredFounder(genome, prob.snp = 1)[2,]

# Recombination event
set.seed(114)
new.gamete <- hypredRecombine(genome,
                              genomeA = founder1,
                              genomeB = founder2,
                              mutate = TRUE,
                              mutation.rate.snp = (2.5 * 10^-5),
                              mutation.rate.qtl = (2.5 * 10^-5),
                              block = FALSE)

# Check the first 20 loci
new.gamete[1:20]

# Let's create a F_2 population
# size: N = 250 > this means that 500 gametes need to be produced
## (each i indidual is the result of a combination of 2 gametes, so 500 gametes make 250 individuals)
N <- 250 # number of individuals
N.g <- N * 2 # number of gametes

# Call on hypredRecombine to generate gametes
F2 <- t(sapply(X = 1:N.g, FUN = function(i) {
  hypredRecombine(genome,
                  genomeA = founder1,
                  genomeB = founder2,
                  mutate = TRUE,
                  mutation.rate.snp = (2.5 * 10^-5),
                  mutation.rate.qtl = (2.5 * 10^-5),
                  block = FALSE)
}))

# This matrix represents the chromosomal sets of each individual
# Rows 1 and 2 are for individual 1, rows 3 and 4 for individual 2, etc

# Random mating
# This is done through nested loops
G <- 50 # number of generations
# Store the F2 gamete information in a matrix
random.mate.pop <- F2
random.mate.pop.tmp <- matrix(nrow = N.g, ncol = 200)

for(generation in 1:G) { ## Loop over generations
  gamete.index.1 <- 1 # Indexing
  gamete.index.2 <- 2
  
  for (indiv in 1:N) { ## Loop over individuals
    random.mate.pop.tmp[gamete.index.1,] <- # For gamete 1
      hypredRecombine(genome,
                      genomeA = random.mate.pop[gamete.index.1,],
                      genomeB = random.mate.pop[gamete.index.2,],
                      mutate = TRUE,
                      mutation.rate.snp = (2.5 * 10^-5),
                      mutation.rate.qtl = (2.5 * 10^-5),
                      block = FALSE)
    
    random.mate.pop.tmp[gamete.index.2,] <- # For gamete 2
      hypredRecombine(genome,
                      genomeA = random.mate.pop[gamete.index.1,],
                      genomeB = random.mate.pop[gamete.index.2,],
                      mutate = TRUE,
                      mutation.rate.snp = (2.5 * 10^-5),
                      mutation.rate.qtl = (2.5 * 10^-5),
                      block = FALSE)
    
    # Increment to the next individual
    gamete.index.1 <- gamete.index.1 + 2
    gamete.index.2 <- gamete.index.2 + 2
    
  } # Close individual loop
  
  # Permutate
  random.mate.pop <- random.mate.pop.tmp[sample(1:N.g),]
  
} # Close generation loop


# Selection
## Obtain genotype values
g.values <- hypredTruePerformance(genome,
                                  random.mate.pop,
                                  DH = FALSE)

# Phenotype values
var.env <- var(g.values)
phen.values <- g.values + rnorm(n = N, 0, sqrt(var.env))

# Select the individual with the highest phenotype
index.individual <- which.max(phen.values)
index.row1 <- (index.individual * 2) - 1
index.row2 <- (index.individual * 2)

selected <- random.mate.pop[c(index.row1, index.row2),]


# Doubled haploid
# This is produced similarly to the F2, except 1 gamete is produced per inidividual
DH <- t(sapply(X = 1:N, FUN = function(i) {
  hypredRecombine(genome,
                  genomeA = selected[1,],
                  genomeB = selected[2,],
                  mutate = TRUE,
                  mutation.rate.snp = (2.5 * 10^-5),
                  mutation.rate.qtl = (2.5 * 10^-5),
                  block = FALSE)
}))

DH.g.values <- hypredTruePerformance(genome,
                                     DH,
                                     DH = T)

selected[, qtl.ids]
selected[, qtl.dom.ids]

# Effect a = 1 times the number of 1 alleles
sum(selected[, qtl.ids] == 1) +
  ## 0.5 if the dominance QTL has genotype 10 or 01
  (sum(selected[, qtl.dom.ids[1]]) == 1) * 0.5 +
  (sum(selected[, qtl.dom.ids[2]]) == 1) * 0.5

# Generate phenotypic values
DH.phen.values <- DH.g.values + rnorm(250, 0, sqrt(var.env))


## Create design matricies
design.DH <- hypredCode(genome,
                        genotypes = DH,
                        DH = T,
                        type = "-101")

design.DH[1, 1:5]
design.DH[1, 17:21]


## Genomic selection (not part of vignette)
corel <- numeric()

for (i in 1:500) {
  train <- sort(sample(1:N, 200))
  validate <- setdiff(1:N, train)
  # Phenos
  y <- DH.phen.values
  y.pred <- as.matrix(y[train,])
  y.validate <- as.matrix(y[validate,])
  
  K <- A.mat(design.DH, min.MAF = 0, max.missing = 1)
  X <- matrix(1, nrow = nrow(y.pred))
  Z <- diag(nrow(K))
  Z <- as.matrix(Z[-validate,])
  
  emm.out <- emmreml(y = y.pred, X = X, Z = Z, K = K, PEVuhat = T)
  
  GEBV.validate <- emm.out$uhat[validate,]
  corel[i] <- cor(GEBV.validate, y.validate)
}



#########
## Testing with barley-relevant parameters
# Initial parameters
n.chr = 7
n.chr.snps = 300
n.QTL = 50


# Create the genome
hv.genome <- hypredGenome(num.chr = n.chr, len.chr = rep(1.5, n.chr), num.snp.chr = n.chr.snps)

# Add QTL
qtl.ids <- seq(from = 60, to = (n.chr * n.chr.snps), by = 60)

hv.genome <- hypredNewQTL(hv.genome,
                          new.id.add = qtl.ids,
                          new.eff.add = rep(1, length(qtl.ids))
                          )

# Create the founder lines
founder1 <- hypredFounder(hv.genome, prob.snp = 1)[1,]
founder2 <- hypredFounder(hv.genome, prob.snp = 1)[2,]

# Create an F_2 population
N = 1000 # Number of individuals
N.gamete = N * 2

F2 <- t(sapply(X = 1:N.gamete, FUN = function(i) {
  hypredRecombine(hv.genome,
                  genomeA = founder1,
                  genomeB = founder2,
                  mutate = TRUE,
                  mutation.rate.snp = (2.5 * 10^-5),
                  mutation.rate.qtl = (2.5 * 10^-5),
                  block = FALSE)
}))

# Inbreed
n.S = 10 # Number of inbreeding generations

inbreeding.pop <- F2
inbreeding.pop.tmp <- inbreeding.pop


    
    

for(generation in 1:G) { ## Loop over generations
  gamete.index.1 <- 1 # Indexing
  gamete.index.2 <- 2
  
  for (indiv in 1:N) { ## Loop over individuals
    random.mate.pop.tmp[gamete.index.1,] <- # For gamete 1
      hypredRecombine(genome,
                      genomeA = random.mate.pop[gamete.index.1,],
                      genomeB = random.mate.pop[gamete.index.2,],
                      mutate = TRUE,
                      mutation.rate.snp = (2.5 * 10^-5),
                      mutation.rate.qtl = (2.5 * 10^-5),
                      block = FALSE)
    
    random.mate.pop.tmp[gamete.index.2,] <- # For gamete 2
      hypredRecombine(genome,
                      genomeA = random.mate.pop[gamete.index.1,],
                      genomeB = random.mate.pop[gamete.index.2,],
                      mutate = TRUE,
                      mutation.rate.snp = (2.5 * 10^-5),
                      mutation.rate.qtl = (2.5 * 10^-5),
                      block = FALSE)
    
    # Increment to the next individual
    gamete.index.1 <- gamete.index.1 + 2
    gamete.index.2 <- gamete.index.2 + 2
    
  } # Close individual loop
  
  # Permutate
  random.mate.pop <- random.mate.pop.tmp[sample(1:N.g),]
  
} # Close generation loop

    
    
    
    
    
    
    