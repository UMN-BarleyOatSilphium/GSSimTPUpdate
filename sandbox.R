## Sandbox

# Let's try to speed this up using rQTL and simcross
library(qtl)
library(simcross)
library(GSsim.TPUpdate)
library(dplyr)
library(stringr)

data("CAP.markers")
data("CAP.haploids")
data("CAP.genotypes")
data("CAP.genos")
data("CAP.marker.map")





CAP.lines <- colnames(CAP.genotypes)[-c(1:4)]

# Crossing block
crossing.block <- make.crossing.block(parent1.lines = CAP.lines[1:100], 
                                      parent2.lines = CAP.lines[1:100], n.crosses = 50,
                                      method = "random", use.parents.once = T)

# Population size per cross
n.ind = 20
# Number of crosses
n.crosses = 10


# Sample founder genotypes
founder1 <- CAP.haploids[c(1,2),] %>%
  apply(MARGIN = 2, FUN = sum) + 1
founder2 <- CAP.haploids[c(1471, 1472),] %>%
  apply(MARGIN = 2, FUN = sum) + 1

founders <- rbind(founder1, founder2) %>% t()

map <- CAP.marker.map

n.self.gen <- 2

m = 0

p = 0

parent.genos <- CAP.genos

ped <- make.pedigree(n.self.gen = n.self.gen, n.ind = n.ind)

# IDs to extract
prog.ids <- which(ped$gen == (n.self.gen + 1))


# Iterate over the crossing block
prog.genos <- apply(X = crossing.block, MARGIN = 1, FUN = function(cross) {
  
  # Extract the parent names
  p1 <- cross[1] %>% as.matrix() %>% as.character()
  p2 <- cross[2] %>% as.matrix() %>% as.character()
  
  # Extract the parent genos
  cross.genos <- parent.genos[c(p1, p2),] %>%
    t()
  
  # Simulate cross-over data
  xo.data <- sim_from_pedigree_allchr(pedigree = ped, map = map, m = 0)
  
  # Generate genotypes
  convert2geno_allchr(xodat = xo.data, map = map, 
                      founder_geno = cross.genos, return.matrix = T,
                      id = prog.ids) 
})




## Make the population
pop1 <- make.population1(map = map, crossing.block = crossing.block, 
                         parent.genos = CAP.genos, n.ind = n.ind, 
                         n.self.gen = n.self.gen, m = 0)
                         




xo.data <- sim_from_pedigree_allchr(pedigree = ped1, map = map, m = 0)







# Extract info for 1st ped
xo.data1 <- lapply(xo.data, FUN = function(chr) {
  to.keep <- startsWith(names(chr), "1")
  chr[to.keep] })

# Generate genos
genos <- convert2geno_allchr(xodat = xo.data1, map = map, return.matrix = T)

# Create a population
population <- make.population(map = CAP.map, crossing.block = crossing.block,
                              parent.genos = CAP.genos, n.ind = 20, n.gen = 2, 
                              cycle.number = 1)

# Make a family of inbred lines
progeny <- make.family2(map = CAP.map, parent1.genome = founder1, 
                        parent2.genome = founder2, n.gen = 2, n.ind = 30, m = 0,
                        cycle.number = 1, family.number = 1)







check_pedigree(ped.all, ignore_sex = T)

xodat <- sim_from_pedigree_allchr(pedigree = ped.all, map = CAP.map, m = 0)

prog.genos <- convert2geno_allchr(xodat = xodat, map = CAP.map, founder_geno = founders, 
                                  id = 23:42, return.matrix = T)


prog.genos <- replicate(n = n.ind, expr = {

  xodat.i <- sim_from_pedigree_allchr(pedigree = ped, map = CAP.map, m = 0)
  
  # Get genotypic data
  convert2geno_allchr(xodat = xodat.i, map = CAP.map, return.matrix = T, id = c(4,5),
                      founder_geno = founders)
  
  }) %>%
  t()
  
  # bind the rows
  do.call("rbind", .) })




xodat.i <- sim_from_pedigree_allchr(pedigree = ped, map = CAP.map, m = 0)

prog.genos <- convert2geno_allchr(xodat = xodat.i, map = CAP.map, 
                                  founder_geno = founders, return.matrix = T, 
                                  id = 4)




table(prog.genos) / sum(table(prog.genos))


# Iterate over the number of crosses
pop.genos <- lapply(X = seq(n.crosses), FUN = function(n) {
  
  # To create a population, we create a list of xodat for each individual
  lapply(X = seq(n.ind), FUN = function(i) {
    xodat.i <- sim_from_pedigree_allchr(pedigree = ped, map = CAP.map, m = 0)
    
    # Get genotypic data
    convert2geno_allchr(xodat = xodat.i, map = CAP.map, return.matrix = T, id = 5)
    
  })
    
    # bind the rows
    do.call("rbind", .) })
  
  # Convert to genotype data
  


