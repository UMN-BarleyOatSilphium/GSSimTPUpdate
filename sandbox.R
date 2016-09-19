## Sandbox

# Let's try to speed this up using rQTL and simcross
library(qtl)
library(simcross)
library(GSsim.TPUpdate)
library(dplyr)

data("CAP.markers")
data("CAP.haploids")
data("CAP.genotypes")


CAP.lines <- colnames(CAP.genotypes)[-c(1:4)]

# Crossing block
crossing.block <- make.crossing.block(parent1.lines = CAP.lines[1:40], 
                                      parent2.lines = CAP.lines[1:40], n.crosses = 20,
                                      method = "random", use.parents.once = T)

# Create a RIL pedigree matrix
ped <- sim_ril_pedigree(ngen = 3, selfing = T, parents = c(1,2))

ped1 <- data.frame(id = c(1,2), mom = c(0,0), dad = c(0,0), sex = c(0,1), gen = c(0,0))

# Population size per cross
n.ind = 20
# Number of crosses
n.crosses = 30

ids = paste(1, seq(n.ind), sep = "") %>% as.integer()

f1.ped <- data.frame(id = ids, mom = rep(1,n.ind), dad = rep(2, n.ind), sex = 0, gen = 1)

ids = paste(2, seq(n.ind), sep = "") %>% as.integer()

f2.ped <- data.frame(id = ids, mom = f1.ped$id, dad = f1.ped$id, sex = 0, gen = 2)

ped.all <- rbind(ped1, f1.ped, f2.ped)

# Genetic map
CAP.map.split <- split(x = CAP.markers, f = CAP.markers$chrom)
CAP.map <- lapply(X = CAP.map.split, FUN = function(chr) {
  chr.map <- chr$pos * 100
  names(chr.map) <- chr$rs
  return(chr.map) })

# Sample founder genotypes
founder1 <- CAP.haploids[c(1,2),] %>%
  apply(MARGIN = 2, FUN = sum) + 1
founder2 <- CAP.haploids[c(1471, 1472),] %>%
  apply(MARGIN = 2, FUN = sum) + 1

founders <- rbind(founder1, founder2) %>% t()

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
  


