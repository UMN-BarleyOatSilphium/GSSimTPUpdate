## Plot the parsed results of the simulation

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


## Define a function for summarization and calculation of mean and confidence
## interval
# df - a data.frame of the update method, cycle number, and parameter value.
sim.summarize <- function(df) {
  
  # Remove duplicated change-iter-cycle pairs
  df1 <- df %>% 
    group_by(change, iter, cycle) %>% 
    filter(row_number() == 1)
  
  df2 <- df1 %>%
    group_by(change, cycle) %>%
    summarize(mean = mean(value, na.rm = T), 
              sd = sd(value, na.rm = T),
              n = n()) %>%
    mutate(se = sd / sqrt(n),
           ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)
  
  # Add in an offset to the cycles
  df2$cycle.offset <- df2$cycle + (as.numeric(as.factor(df2$change)) - 1) * 0.1
  
  return(df2)
}




## Define a function to plot the summary data.frame using ggplot2

sim.ggplot <- function(df.summary, main, ylab, xlab = "Breeding Cycle", ylim, 
                       xlim = c(1, n.cycles + 1), tp.change.factors) {
  
  gp <- ggplot(data = df.summary, aes(x = cycle.offset, y = mean, col = change, 
                                      shape = change)) +
    geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    ggtitle(main) + 
    ylab(ylab) +
    xlab(xlab) +
    ylim(ylim) +
    xlim(xlim) +
    scale_color_discrete(name = "Update Method",
                         labels = as.character(tp.change.factors)) +
    scale_shape_discrete(name = "Update Method",
                         labels = as.character(tp.change.factors))
  
  print(gp)
  
}
  


### Change in Genetic Variance
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$candidate.gen.var) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, main = "Genetic Variance", ylab = "Genetic Variance",
           ylim = c(0,8), tp.change.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_gen_var.jpg", sep = "")),
       height = 5, width = 6.5)


### Change in mean genotypic value of the selection candidates

df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$candidate.gen.val) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, main = "Genotypic Value of the Selection Candidates", 
           ylab = "Genotypic Value", ylim = c(0,25), tp.change.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_geno_value.jpg", sep = "")),
       height = 5, width = 6.5)

## Alternatively, the true response to selection

df1 <- df %>% 
  group_by(change, iter, cycle) %>% 
  filter(row_number() == 1) %>%
  group_by(change, iter) %>%
  mutate(value = c(NA, diff(value))) %>% 
  na.omit() %>%
  sim.summarize()

sim.ggplot(df.summary = df1, main = "Response to Selection Across Cycles", 
           ylab = "Response to Selection", ylim = c(-1,5), 
           tp.change.factors = tp.change.factors)


### Change in the prediction accuracy over cycles

# Plot
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$validation.results) %>%
  bind_rows()

df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, main = "Realized Prediction Accuracy", 
           ylab = "Prediction Accuracy (r)", ylim = c(0,0.75), 
           tp.change.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_prediction_accuracy.jpg", sep = "")),
       height = 5, width = 6.5)



## Potential Explanatory Variables


### Change in QTL-marker LD over cycles

# Each QTL's max LD
# Plot
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$qtl.marker.LD) %>%
  bind_rows()

df1 <- sim.summarize(df %>% filter(extra == "mean_max_genome"))
  

sim.ggplot(df.summary = df1, main = "Mean LD Between QTL and Marker in Highest LD", 
           ylab = "Linkage Disequilibrum (r)", ylim = c(0.6, 0.9), 
           tp.change.factors = tp.change.factors)


ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_mean_max_ld.jpg", sep = "")),
       height = 5, width = 6.5)



## Persistence of phase - use the same df

df1 <- sim.summarize(df %>% filter(extra == "persistence"))


sim.ggplot(df.summary = df1, 
           main = "Persistence of LD Phase Between Training Population\nand Selection Candidates", 
           ylab = "Persistence of Phase (cor of r)", ylim = c(0,0.65), 
           tp.change.factors = tp.change.factors)


ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_persistence_of_phase.jpg", sep = "")),
       height = 5, width = 6.5)

### Relationship
## TP - SC

# Plot
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$tp.sc.relationship) %>%
  bind_rows()

df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, 
           main = "Additive Genomic Relationship Between Training Population\nand Selection Candidates", 
           ylab = "Additive Genomic Relationship\n(Scaled to Base Population)", 
           ylim = c(-0.1,1.25), tp.change.factors = tp.change.factors)


ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_average_relationship.jpg", sep = "")),
       height = 5, width = 6.5)


### Change in inbreeding coefficient
## Selection candidates

df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$sc.inbreeding) %>%
  bind_rows()

df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, 
           main = "Average Inbreeding Coefficient\nAmong Selection Candidates", 
           ylab = "", 
           ylim = c(0, 2.5), tp.change.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_sc_inbreeding.jpg", sep = "")),
       height = 5, width = 6.5)


## Rate of inbreeding
df1 <- df %>% 
  group_by(change, iter, cycle) %>% 
  filter(row_number() == 1) %>%
  group_by(change, iter) %>%
  mutate(value = c(NA, diff(value))) %>% 
  na.omit() %>%
  sim.summarize()

sim.ggplot(df.summary = df1, 
           main = "Rate of Inbreeding\nAmong Selection Candidates", 
           ylab = "Rate of Inbreeding (delta F)", 
           ylim = c(0, 0.3), tp.change.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_sc_inbreeding_rate.jpg", sep = "")),
       height = 5, width = 6.5)



## TP additions

df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$tp.additions.inbreeding) %>%
  bind_rows()

df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, 
           main = "Average Inbreeding Coefficient\nAmong Training Population Additions", 
           ylab = "", 
           ylim = c(0, 2.5), tp.change.factors = tp.change.factors)

## Parents

df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$parent.inbreeding) %>%
  bind_rows()

df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, 
           main = "Average Inbreeding Coefficient\nAmong Breeding Cycle Parents", 
           ylab = "", 
           ylim = c(0, 2.5), tp.change.factors = tp.change.factors)

### Change in the Expected Heterozygosity Among TP Additions

# Plot
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$tp.update.exp.het) %>%
  bind_rows()

df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, 
           main = "Expected Heterozygosity of Training Population Additions", 
           ylab = "Expected Heterozygosity", ylim = c(0,0.25), 
           tp.change.factors = tp.change.factors)


ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_exp_het.jpg", sep = "")),
       height = 5, width = 6.5)



## Changes in QTL allele frequency

# Proportion of fixed QTL in the SC
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$candidate.prop.fixed) %>%
  bind_rows()

df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, 
           main = "Proportion of QTL Fixed for an Allele in the Selection Candidates", 
           ylab = "Proportion", ylim = c(0,1), 
           tp.change.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_fixed_qtl.jpg", sep = "")),
       height = 5, width = 6.5)


# # Save data lists to a file
# # Gather the lists
# data.lists <- neyhart::find('\\.list$', class = "list")
# 
# # Create a new filename
# save.file <- filename %>%
#   str_replace(pattern = "simulation_results", replacement = "plot_data")
# 
# # Save the data
# save(list = data.lists, file = save.file)
