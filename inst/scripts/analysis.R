## Script to perform analysis of simulation results

# Project: BarleySimGS-TPUpdate
# Author: Jeff Neyhart
# Date: 16 May 2016

# Load packages
library(dplyr)
library(stringr)

# Designate the directory with the files
results.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/BarleySimGS-TPUpdate/Results/Base Experiment/"

# Load data from the allele frequency experiment
all.files <- list.files(results.dir, full.names = T)
filename <- all.files[1]

load(filename)

# Experimental parameters
pop.makeup <- str_extract(string = filename, pattern = 'popmakeup-[A-Za-z]{2,5}') %>%
  str_extract(pattern = '[A-Za-z]{2,5}$')

tp.formation <- str_extract(string = filename, pattern = 'tpformation-[a-z]*') %>% 
  str_extract(pattern = '[a-z]*$')

n.cycles = collective.abbreviated.results[[1]][[1]][[1]][[1]] %>%
  length()
n.reps = lapply(X = collective.abbreviated.results[[1]][[1]], FUN = length) %>%
  unlist() %>%
  sum()


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




