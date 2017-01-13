## Gather and Summarize Parsed Results for Plotting

# Load packages
library(tidyverse, quietly = T)
library(stringr, quietly = T)
library(GSSimTPUpdate)

# Set directory and grab the files
results.dir <- "/panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Simulations/GSSimTPUpdate/inst/output"

# Load data
all.files <- list.files(results.dir, full.names = T, pattern = "collective")

# Create a list to store multiple collective data list
total.collective.data <- list()

# Iterate over files and assign data to designators
for (file in all.files) {
  # Grab some metadata
  pop.type <- file %>% basename() %>% str_extract('MNxND|MN|ND')
  # TP formation
  tp.formation <- file %>% basename() %>% str_extract('cumulative|window')
  # Heritability
  h2 <- file %>% basename %>% str_extract('05|02')
  
  load(file)
  
  # Assign the data to a new object
  assign.name <- str_c(pop.type, tp.formation, h2, sep = "_")
  total.collective.data[[assign.name]] <- collective.abbreviated.results
  
}

# Get rid of one list
rm(collective.abbreviated.results)

tp.change.factors <- c(best = "Top", CDmean = "CDmean", nochange = "No Change", 
                       PEVmean = "PEVmean", random = "Random", worst = "Bottom") %>%
  as.factor()

# Number of cycles
n.cycles <- total.collective.data[[1]][[1]][[1]]$cycle %>% unique() %>% length()


# Names of the larger list
total.names <- names(total.collective.data)

# Create a list to store the data.frames
df.list <- list()

# Prediction accuracy
df.list[["accuracy"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$validation.results) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()

# Genetic variance

df.list[["genvar"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.gen.var) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()

# Genotypic value

df.list[["genval"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.gen.val) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()

# QTL - marker LD

df.list[["qtl_marker_LD"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$qtl.marker.LD) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()

# Average relationship

# Plot
df.list[["rel"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$tp.sc.relationship) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()


# Inbreeding

df.list[["inbreeding"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.inbreeding) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()

# Fixed QTL

df.list[["fixedqtl"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.allele.freq) %>%
    bind_rows() %>%
    filter(extra1 == "qtl") %>% 
    group_by(change, iter, cycle) %>% 
    summarize(value = sum(value == 0 | value == 1)) %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()

# Fixed markers

df.list[["fixedmar"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.allele.freq) %>%
    bind_rows() %>%
    filter(extra1 == "markers") %>% 
    group_by(change, iter, cycle) %>% 
    summarize(value = sum(value == 0 | value == 1)) %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()

# QTL effects

# test <- lapply(X = total.names, FUN = function(coll.name)
#   lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) {
#     # # QTL effects
#     # eff <- tpc$qtl.effects %>%
#     #   rename(effect = value)
#     # Gather the allele frequencies
#     freq <- tpc$sc.allele.freq %>%
#       filter(extra1 == "qtl")
#     # Combine
#     # full_join(freq, eff)
#     })) %>%
#     bind_rows() %>%
#     mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
#            heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
#   bind_rows()
# 
# 
# 
# 
# 
#     # Gather the effects
#     eff <- tpc$qtl.effects
#     # Add the effects to the frequencies
#     freq$effect <- eff$value
#     return(freq) }) %>%
#     bind_rows() %>%
#     mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
#            heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
#   bind_rows()

# Expected het

df.list[["exphet"]] <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$tp.update.exp.het) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()


# Return the data
save("df.list", file = "data_to_plot.RData")