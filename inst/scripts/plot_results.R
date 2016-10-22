## Plot the parsed results of the simulation

# Load packages
library(tidyverse, quietly = T)
library(stringr, quietly = T)
library(GSsim.TPUpdate)


# Set directory and grab the files
results.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/GSsim.TPUpdate/Results"

# File path for figures
figures.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/GSsim.TPUpdate/figures/"

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


### Figure 3 is the change in prediction accuracy

df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$validation.results) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(coll.name, 'window|cumulative') %>% str_to_title(),
           heritability = str_extract(coll.name, '05|02') %>% str_replace("0", "0.")) ) %>%
  bind_rows()

df1.acc <- sim.summarize(df) %>%
  mutate(variable = "Prediction Accuracy")

sim.ggplot(df.summary = df1.acc, main = "Realized Prediction Accuracy", 
           ylab = expression(Prediction~Accuracy~(italic(r[MG]))), 
           col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_pred_acc_combined.jpg")),
       height = 5, width = 9)


### Figure 4 is the change in genetic variance and genetic value

### Change in Genetic Variance
df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.gen.var) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1.genvar <- sim.summarize(df) %>%
  mutate(variable = "Genetic Variance")

### Change in mean genotypic value of the selection candidates

df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.gen.val) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1.genval <- sim.summarize(df) %>%
  mutate(variable = "Genotypic Value")


### Create a figure with both plots
df1 <- bind_rows(df1.genval, df1.genvar)

# Plot
sim.ggplot(df.summary = df1, main = "", ylab = "", col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_gen_val_gen_var_combined.jpg")),
       height = 8, width = 9)



## Alternatively, the true response to selection

df1 <- df %>% 
  group_by(exp_name, change, iter, cycle) %>% 
  filter(row_number() == 1) %>%
  group_by(exp_name, change, iter) %>%
  mutate(value = c(NA, diff(value))) %>% 
  na.omit() %>%
  ungroup() %>%
  sim.summarize() %>%
  # Add the variable designator
  mutate(variable = "Response to Selection")

sim.ggplot(df.summary = df1, main = "Response to Selection Across Cycles", 
           ylab = "Response to Selection", 
           col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_resp_sel.jpg")),
       height = 5, width = 9)



### Figure 5 will have the genomic relationship, persistence of phase, inbreeding,
### and number of QTL fixed for an allele

### Change in QTL-marker LD over cycles

# Each QTL's max LD

## In the TP
df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$qtl.marker.LD) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()

df1.mean.max <- sim.summarize(df %>% filter(extra1 == "tp_mean_max_genoms")) %>%
  mutate(variable = "Mean Max LD in Training Population")
  

gp.max.LD.tp <- sim.ggplot(df.summary = df1.mean.max, main = "Mean LD Between QTL and Marker in Highest LD", 
           ylab = expression(Linkage~Disequilibrium~(r^2)), text.y.scaling = 1.01,
           col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_tp_mean_max_LD.jpg")),
       height = 5, width = 9)


# In the candidates
df1.mean.max <- sim.summarize(df %>% filter(extra1 == "sc_mean_max_genome")) %>%
  mutate(variable = "Mean Max LD in Selection Candidates")


gp.max.LD.sc <- sim.ggplot(df.summary = df1.mean.max, 
                           main = "Mean LD Between QTL and Marker in Highest LD", 
                        ylab = expression(Linkage~Disequilibrium~(r^2)), text.y.scaling = 1.01,
                        col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_sc_mean_max_LD.jpg")),
       height = 5, width = 9)




## Persistence of phase - use the same df

df1.persistence <- sim.summarize(df %>% filter(extra1 == "persistence")) %>%
  mutate(variable = "Persistence of LD Phase")


### Relationship
## TP - SC

# Plot
df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$tp.sc.relationship) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()

df1.relatioship <- sim.summarize(df) %>%
  mutate(variable = "Average Relationship")


# Inbreeding
df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.inbreeding) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()

df1.inbreeding <- sim.summarize(df) %>%
  mutate(variable = "Inbreeding")

# Fixed QTL
df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.allele.freq) %>%
    bind_rows() %>%
    filter(extra1 == "qtl") %>% 
    group_by(change, iter, cycle) %>% 
    summarize(value = sum(value == 0 | value == 1)) %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()

df1.qtlfreq <- sim.summarize(df %>% ungroup()) %>%
  mutate(variable = "Number of Fixed QTL")




df1 <- bind_rows(df1.persistence, df1.relatioship, df1.inbreeding, df1.qtlfreq)

gp.combined <- sim.ggplot(df.summary = df1, 
                          main = "", 
                          ylab = "", 
                          col.factors = tp.change.factors, text.y.scaling = 1.01)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_explan_vars.jpg")),
       plot = gp.combined, height = 12, width = 9)


### Change in inbreeding coefficient
## Selection candidates

df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.inbreeding) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()

df1.inbreeding <- sim.summarize(df) %>%
  mutate(variable = "Inbreeding")

sim.ggplot(df.summary = df1.inbreeding, 
           main = "Average Inbreeding Coefficient\nAmong Selection Candidates", 
           ylab = "Inbreeding Coefficient", col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_inbreeding_combined.jpg")),
       height = 5, width = 9)


## Rate of inbreeding
df1 <- df %>% 
  group_by(exp_name, change, iter, cycle) %>% 
  filter(row_number() == 1) %>%
  group_by(exp_name, change, iter) %>%
  mutate(value = c(NA, diff(value))) %>% 
  na.omit() %>%
  ungroup() %>%
  sim.summarize() %>%
  mutate(variable = "Rate of Inbreeding")

sim.ggplot(df.summary = df1, 
           main = "Rate of Inbreeding\nAmong Selection Candidates", 
           ylab = "Rate of Inbreeding (delta F)", col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_inbreeding_rate_combined.jpg")),
       height = 5, width = 9)


### Change in the Expected Heterozygosity Among TP Additions

# Plot
df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$tp.update.exp.het) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()


df1.exphet <- sim.summarize(df) %>%
  mutate(variable = "Expected Heterozygosity")

sim.ggplot(df.summary = df1.exphet, 
           main = "Expected Heterozygosity of Training Population Additions", 
           ylab = "Expected Heterozygosity", 
           col.factors = tp.change.factors)


ggsave(filename = file.path(figures.dir, str_c(pop.type, "_exp_het_combined.jpg")),
       height = 5, width = 9)



## Changes in QTL allele frequency

# Proportion of fixed QTL in the SC
df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.allele.freq) %>%
    bind_rows() %>%
    filter(extra1 == "qtl") %>% 
    group_by(change, iter, cycle) %>% 
    summarize(value = sum(value == 0 | value == 1)) %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%
  bind_rows()

df1.qtlfreq <- sim.summarize(df %>% ungroup())

sim.ggplot(df.summary = df1, 
           main = "Number of QTL Fixed for an Allele in the Selection Candidates", 
           ylab = "Number of QTL", 
           col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_prop_fixed_combined.jpg")),
       height = 5, width = 9)


## Combine QTL with their effect in each iteration
## Find the proportion of QTL that are fixed for the
## favorable allele
df <- lapply(X = total.names, FUN = function(coll.name) 
  lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) {
    # Gather the allele frequencies
    freq <- tpc$sc.allele.freq %>%
      filter(extra1 == "qtl")
    # Gather the effects
    eff <- tpc$qtl.effects
    # Add the effects to the frequencies
    freq$effect <- eff$value
    return(freq) }) %>%
    bind_rows() %>%
    mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
             str_to_title()) ) %>%

  bind_rows() %>%
  # Find QTL fixed for favorable allele
  group_by(exp_name, change, iter, cycle) %>% 
  filter((value == 0 & effect < 0) | (value == 1 & effect > 0) ) %>%
  summarize(value = n())
  
df1 <- sim.summarize(df %>% ungroup())
    
sim.ggplot(df.summary = df1, 
           main = "Number of QTL Fixed for the Favorable Allele in the Selection Candidates", 
           ylab = "Number of QTL", 
           col.factors = tp.change.factors)         
         
         

