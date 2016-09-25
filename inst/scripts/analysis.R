## Analysis of data

# Load packages
library(tidyverse, quietly = T)
library(stringr, quietly = T)
library(GSsim.TPUpdate, quietly = T, verbose = F)


# Set directory and grab the files
results.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/GSsim.TPUpdate/Results/Base Experiment/"

# File path for figures
figures.dir <- file.path(results.dir, "Figures")

# Load data from the allele frequency experiment
all.files <- list.files(results.dir, full.names = T) %>%
  str_subset(pattern = "simulation_results") %>%
  str_subset(pattern = "collective")

# Just MNxND files
filename <- all.files %>%
  str_subset(pattern = "MNxND") %>%
  .[1]

load(filename)

### Initial experimen settings
# Experimental parameters
pop.makeup <- str_extract(string = filename, pattern = 'results_[A-Za-z]{2,5}') %>%
  str_extract(pattern = '[A-Za-z]{2,5}$')

tp.formation <- 
  str_extract(string = filename, pattern = paste0(c("cumulative", "window"), collapse = "|"))

n.cycles <- collective.abbreviated.results[[1]][[1]] %>% 
  select(cycle) %>% 
  unique() %>% 
  nrow()

# The different updating methods
collective.abbreviated.results %>%
  names()

tp.change.factors <- c(best = "Best", CDmean = "CDmean", nochange = "No Change", 
                       PEVmean = "PEVmean", random = "Random", worst = "Worst") %>%
  as.factor()



## Try to model prediction accuracy as a function of genetic variance and 
## persistence of LD phase

# Gather data on these three parameters
multi.df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
  list(tpc$validation.results, 
       tpc$candidate.gen.var, 
       filter(tpc$qtl.marker.LD, extra == "persistence")) %>% 
    reduce(full_join, by = c("change", "iter", "cycle")) %>%
    rename(accuracy = value.x, V_g = value.y, persistence = value) %>%
    select(-extra) )

# Model

  


# Find markers in each replication of the experiment that have gone from a
## low MAF to an intermediate MAF
drifted.marker.list <- lapply(collective.abbreviated.results, function(tpc)
  unlist(lapply(X = tpc$allele.freq.list, FUN = function(set)
    sapply(X = set, FUN = function(rep) {
      # Find the markers in cycle 1 that are at a low frequency
      maf.c1 <- sapply(rep$cycle1, function(freq) min(freq, 1 - freq))
      mar.c1 <- names(which(maf.c1 > 0 & maf.c1 < 0.05))
      # Find the markers in cycle 15 that are at a higher frequency
      maf.c15 <- sapply(rep$cycle15, function(freq) min(freq, 1 - freq))
      mar.c15 <- names(which(maf.c15 > 0.2 & maf.c15 < 0.8))
      # Find the markers in common
      intersect(mar.c1, mar.c15) })), 
    recursive = F) )




