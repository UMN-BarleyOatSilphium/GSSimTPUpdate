## Sandbox

# Let's try to speed this up using rQTL and simcross
library(qtl)
library(simcross)
library(GSsim.TPUpdate)
library(dplyr)
library(mpMap2)

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

map <- CAP.marker.map

n.self.gen <- 2

m = 0

p = 0

parent.genos <- CAP.genos

ped1 <- ped2 <- make.pedigree(n.self.gen = n.self.gen, n.ind = n.ind)

ped1 <- ped1 %>% 
  mutate(id = str_c(1,id) %>% as.numeric()) %>%
  mutate(mom = str_c(1,mom) %>% as.numeric()) %>%
  mutate(dad = str_c(1,dad) %>% as.numeric())
ped1[1:2,2:3] <- 0

ped2 <- ped2 %>% 
  mutate(id = str_c(2,id) %>% as.numeric()) %>%
  mutate(mom = str_c(2,mom) %>% as.numeric()) %>%
  mutate(dad = str_c(2,dad) %>% as.numeric())
ped2[1:2,2:3] <- 0

ped <- rbind(ped1, ped2) %>%
  arrange(id)

xo.data <- sim_from_pedigree_allchr(pedigree = ped, map = map, m = 0)

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
  


