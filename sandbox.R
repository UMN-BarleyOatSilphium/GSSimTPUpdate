## Sandbox

# Let's try to speed this up using rQTL and simcross
library(qtl)
library(simcross)
library(GSsim.TPUpdate)
library(dplyr)

# Create a RIL pedigree matrix
ped <- sim_ril_pedigree(ngen = 2, selfing = T, parents = c(1,2))

# Population size per cross
n.ind = 20
# Number of crosses
n.crosses = 50

# Genetic map
CAP.map.split <- split(x = CAP.markers, f = CAP.markers$chrom)
CAP.map <- lapply(X = CAP.map.split, FUN = function(chr) {
  chr.map <- chr$pos * 100
  names(chr.map) <- chr$rs
  return(chr.map) })

TP.haploids.i <- CAP.haploids

# Create a crossing block
crossing.block.i <- make.crossing.block(parent1.lines = parent.lines.list$p1, 
                                        parent2.lines = parent.lines.list$p2, 
                                        n.crosses = n.crosses, 
                                        use.parents.once = T)

# Use CAP line data for the founders


# Iterate over the number of crosses
system.time(pop.genos <- lapply(X = seq(n.crosses), FUN = function(n) {
  
  # To create a population, we create a list of xodat for each individual
  lapply(X = seq(n.ind), FUN = function(i) {
    xodat.i <- sim_from_pedigree_allchr(pedigree = ped, map = CAP.map, m = 0)
    
    # Get genotypic data
    convert2geno_allchr(xodat = xodat.i, map = CAP.map, return.matrix = T, id = 5)
    
    # bind the rows
    do.call("rbind", .) }))
  
  # Convert to genotype data
  


